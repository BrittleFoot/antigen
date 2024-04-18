#!python
"""
requirements

```bash
python -m pip install pandas matplotlib biopython dna_features_viewer bcbio-gff
```
"""
import argparse
import io
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import pandas as pd
from BCBio import GFF
from dna_features_viewer import BiopythonTranslator
from matplotlib import pyplot as plt

OANTIGENE_GENES = {
    "wzz",
    "wzx",
    "wzy",
    "wzt",
    "wzm",
    "rfbA",
    "rfbB",
    "rfbC",
    "rfbD",
    "vioA",
    "ugd",
    "galE",
    "manB",
    "manC",
    "gmd",
    "rmd",
    "rfaL",
    "waaL",
    "rmlA",
    "rmlB",
    "rmlC",
    "rmlD",
    "manA",
}


def parse_line(line: str) -> [int, list]:
    row = line.split("\t")

    # contig number
    if len(row) == 1:
        return int(row[0]), None

    return None, row


ID_PATTERN = re.compile(r"ID=(\w+);")


def parse_coordinates(coords_file: io.FileIO) -> pd.DataFrame:
    data = {}
    for record in GFF.parse(coords_file):
        for feature in record.features:
            data[feature.qualifiers["ID"][0]] = record.id

    return data


def fix_numeric(df: pd.DataFrame, key: str):
    df[key] = df[key].apply(pd.to_numeric, errors="coerce").astype("Int64")


def parse_list_of_operons(
    operonfinder_file: io.FileIO,
    coordinates: pd.DataFrame,
) -> pd.DataFrame:
    """Parse the output of operonfinder and coordinates file to a pandas dataframe"""

    operon = 0
    header = None
    data = []

    for line in operonfinder_file:
        line = line.strip()

        if header is None:
            header = line.split("\t")
            header += ["Contig"]
            continue

        operon_update, row = parse_line(line)
        if operon_update is not None:
            operon = operon_update
            continue

        if row is not None:
            # row content:
            # IdGene	Type	COGgene	PosLeft	postRight	Strand	Function
            contig = coordinates.get(row[0])
            if not contig:
                continue

            row = [operon] + row + [contig]

            data.append(row)

    frame = pd.DataFrame(data, columns=header)

    fix_numeric(frame, "Operon")
    fix_numeric(frame, "PosLeft")
    fix_numeric(frame, "postRight")

    return frame


def find_oantigenes(annotation):
    genes = []
    for record in GFF.parse(annotation):
        for feature in record.features:
            q = feature.qualifiers
            gene_names = (q.get("gene", []) + q.get("Gene", [])) or q.get("Name", [])
            if not gene_names:
                continue

            if len(gene_names) > 1:
                raise NotImplementedError(
                    f"More than one gene name: {gene_names} in {feature}"
                )

            gene = gene_names[0]

            # replace to just flat `in` if substring is not needed,
            # and the whole gene in 'gene' qualifier
            good_gene = False
            for oantigene_name in OANTIGENE_GENES:
                if oantigene_name.find(gene) != -1:
                    good_gene = True
                    break

            if not good_gene:
                continue

            q["_record_id"] = record.id

            genes.append(feature)

    return genes


@dataclass
class Operon:
    id: int
    contig: str
    start: int
    end: int
    functions: list


def group_operons(operons_df: pd.DataFrame) -> Dict[int, Operon]:
    groups = dict()

    for i, feature in operons_df.iterrows():
        if feature.Operon not in groups:
            groups[feature.Operon] = Operon(
                id=feature.Operon,
                contig=feature.Contig,
                start=feature.PosLeft,
                end=feature.postRight,
                functions=[feature.Function],
            )
            continue

        op = groups[feature.Operon]

        assert op.contig == feature.Contig, (op, feature)

        op.start = min(op.start, feature.PosLeft)
        op.end = max(op.end, feature.postRight)
        op.functions.append(feature.Function)

    return groups


def find_oantigen_operons(operons: Dict[int, Operon], oantigens: list) -> List[Operon]:
    ids = []
    by_name = []
    for i, operon in operons.items():
        for oantigen in oantigens:
            q = oantigen.qualifiers
            if gff_key(operon.contig) != gff_key(q["_record_id"]):
                continue

            loc = oantigen.location
            if operon.start <= loc.start and loc.end <= operon.end:
                ids.append(operon)

        for func in operon.functions:
            if func.find("O-antigen") != -1:
                by_name.append(operon)

    return ids + by_name


def draw_region_by_coordinates(
    draw_ax,
    operon_number: int,
    start: int,
    end: int,
    record,
    args: "InputArgs",
):
    """
    Function produces schematic plot of DNA region specified with coordinates;
    visualizes only gene feature types
    and labels them by gene or its product.
    """
    labels = [
        "Name",
        "name",
        "product",
        "source",
        "locus_tag",
        "note",
    ]

    if args.label:
        labels = [args.label] + labels

    translator = BiopythonTranslator()
    translator.ignored_features_types = ["region"]
    translator.label_fields = labels
    graphic_record = translator.translate_record(record)

    fig = plt.figure(figsize=(10 * args.scale, 3 * args.scale))
    try:
        operon = graphic_record.crop((start, end))
    except Exception as e:
        print(f"Error on {start, end}, {e}")
        return fig

    ax = fig.add_subplot()

    for a in [draw_ax, ax]:
        if a is None:
            continue

        a.set_title(f"Operon #{operon_number}; Contig: {record.id}")
        operon.plot(ax=a, elevate_outline_annotations=False)

    return fig


NODE_ID_PATTERN = re.compile(r"(\d+)")


def gff_key(string):
    """
    Operon mapper tends to crop and rename contig names with > 21 length,
    preffixes, though, stays the same
    """
    ids = NODE_ID_PATTERN.findall(string)
    if ids:
        return ids[0]

    return string[:21]


def gff_records(annotation: Path):
    return {gff_key(record.id): record for record in GFF.parse(annotation)}


@dataclass
class InputArgs:
    annotation: Path
    operons: Path
    coordinates: Path
    plot: bool
    out_dir: Path | None
    scale: int
    label: str


def parse_args(args: list[str] = None):
    arp = argparse.ArgumentParser("oantigen")

    g = arp.add_argument_group("Annotation file")
    g.add_argument("annotation", type=Path, help="gff annotation file")

    g = arp.add_argument_group(
        "Operon mapper",
        description="The files you will get as a result from the Operon Mapper run",
    )
    g.add_argument("operons", type=Path, help="list_of_operons file")
    g.add_argument("coordinates", type=Path, help="ORFs_coordinates file")

    g = arp.add_argument_group("Output options")
    g.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="use matplotlib to show window with result images (not very friendly if there are more than 6-8 oantigenes)",
    )
    g.add_argument(
        "-o",
        "--out-dir",
        type=Path,
        default=None,
        help="output images directory. by default will write to oantigene_figures_`annotation` directory",
    )
    g.add_argument(
        "-s",
        "--scale",
        type=int,
        default=1,
        help="output images size scale factor. Make it bigger if genes do not fit in the base size",
    )
    g.add_argument(
        "-l",
        "--label",
        default="gene",
        help="which label for the gene to use (from gff qualifiers)",
    )

    parsed = InputArgs(**arp.parse_args(args=args).__dict__)
    if parsed.out_dir is None:
        out_name = "oantigene_figures_" + parsed.annotation.name
        parsed.out_dir = parsed.annotation.parent / out_name

    parsed.out_dir.mkdir(exist_ok=True, parents=True)
    return parsed


def setup_axes(n, real):
    n2 = math.ceil(n / 2)

    axes = [None] * n
    if real:
        _, axes = plt.subplots(
            nrows=n2,
            ncols=2,
            figsize=(14, 2 * n2),
            layout="constrained",
        )
        axes = axes.flatten()

    return axes


def main():
    args = parse_args(
        # "rast/assembly.gff3 rast/list_of_operons_1758265 rast/ORFs_coordinates_1758265".split()
    )

    with args.coordinates.open() as coords_file:
        coordinates = parse_coordinates(coords_file)

    with args.operons.open() as operons_file:
        operons_df = parse_list_of_operons(operons_file, coordinates)

    gene_coordinates = find_oantigenes(args.annotation)
    operons = group_operons(operons_df)
    selected_operons = find_oantigen_operons(operons, gene_coordinates)

    selected_operons = sorted(selected_operons, key=lambda o: o.id)

    records = gff_records(args.annotation)

    axes = setup_axes(len(selected_operons), args.plot)
    for ax, operon in zip(axes, selected_operons):
        contig_key = gff_key(operon.contig)
        record = records[contig_key]

        fig = draw_region_by_coordinates(
            ax,
            operon.id,
            # Do not crop the start of the gene
            operon.start - 1,
            operon.end,
            record,
            args,
        )

        fig.savefig(args.out_dir / f"operon_{operon.id}__{record.id}.png")
        plt.close(fig)

    if args.plot:
        plt.show()


if __name__ == "__main__":
    main()
