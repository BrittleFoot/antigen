#!python
"""
requirements

```bash
python -m pip install pandas matplotlib biopython dna_features_viewer bcbio-gff
```
"""
import argparse
import io
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from BCBio import GFF
from dna_features_viewer import BiopythonTranslator
from matplotlib import pyplot as plt

OANTIGENE_GENES = [
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
]


def parse_line(line: str) -> [int, list]:
    row = line.split("\t")

    # contig number
    if len(row) == 1:
        return int(row[0]), None

    return None, row


ID_PATTERN = re.compile(r"ID=(\w+);")


def parse_coordinates(coords_file: io.FileIO) -> pd.DataFrame:
    data = {}
    for line in coords_file:
        line = line.strip()
        row = line.split("\t")

        contig_name = row[0]
        quals = row[8]

        id = ID_PATTERN.search(quals)
        data[id.group(1)] = contig_name

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


def find_oantigenes_data_from_gff(gff_file):
    """
    Function gets data (start and end coordinates, strand (+,-)) on candidate O-antigen operon genes from GFF annotation
    using list of O-antigen operon genes defined from literature
    :param gff_file: path to GFF annotation file (PGAP or PROKKA)
    :param O_ag_gene_list: list of O-antigen operon genes defined from literature
    :param prokka: True or False; specifies type of GFF annotation
    :return: 2-D array with 4 columns: start coordinate, end coordinate, strand, gene name
    """
    with open(gff_file) as gff_with_fasta:
        i = 0
        start_fasta = 0
        for i, line in enumerate(gff_with_fasta):
            if line.startswith("##FASTA"):
                start_fasta = i
        end_fasta = i

    if start_fasta == 0:
        start_fasta = end_fasta

    gff = pd.read_csv(
        gff_file,
        engine="python",
        sep="\t",
        comment="#",
        skipfooter=end_fasta - start_fasta + 1,
    )

    gff.columns = [i for i in range(1, len(gff.columns) + 1)]
    gff_gene = gff.loc[gff[9].str.contains("Name=")].copy()
    gff_gene.loc[:, 10] = gff_gene[9].apply(lambda x: x.split("Name=")[1].split(";")[0])

    return (
        gff_gene[gff_gene[9].apply(lambda x: any([k in x for k in OANTIGENE_GENES]))]
        .loc[:, [4, 5, 7, 10]]
        .to_numpy()
    )


def find_oantigene_operon_numbers(operons_df, coord_genes_array):
    """
    Function gets candidate O-ag operons (as operon numbers) from Operon-mapper output
    using coordinates array of O-ag genes obtained from GFF annotation
    :param operons_df: pandas DataFrame of Operon-mapper output with the description of predicted operons
    :param coord_genes_array: 2-D array with the following columns: start coordinate, end coordinate, strand, gene name.
    :return: numpy array: operon numbers of candidate O-antigen operons
    """
    return operons_df.loc[
        operons_df.PosLeft.isin(coord_genes_array[:, 0])
        & operons_df.postRight.isin(coord_genes_array[:, 1])
        | operons_df.Function.str.contains("O-antigen")
    ]["Operon"].unique()


def get_operon_boundaries(operons_df, operon_number):
    """
    Function gets boundary coordinates of the DNA region specified as the list of operons.
    :param operons_df: pandas DataFrame of Operon-mapper output with the description of predicted operons
    :param operon_numbers_list: list of operon numbers (int) predicted by Operon-mapper
    :return: tuple consisting of 2 integers:
    (<start_coordinate_of_the_first_gene_from_the_first_operon>, <end_coordinate_of_the_last_gene_from_the_last_operon>)
    """
    operon_boundaries = operons_df[operons_df.Operon.eq(operon_number)]
    operon_boundaries = operon_boundaries.loc[:, ["PosLeft", "postRight", "Contig"]]
    contig = operon_boundaries.Contig.unique()[0]

    s, e = sorted(
        map(int, [operon_boundaries.PosLeft.min(), operon_boundaries.postRight.max()])
    )

    return s, e, contig


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
    :param gff_file: path to GFF annotation file
    :param start_: start coordinate of the desired region
    :param end_: end coordinate of the desired region
    :param prokka: True or False; specifies type of GFF annotation
    :param gff_without_fasta: default None, if Prokka it's used to create new gff
    :return: None
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
    translator.label_fields = labels
    graphic_record = translator.translate_record(record)
    operon = graphic_record.crop((start, end))

    fig = plt.figure(figsize=(10 * args.scale, 3 * args.scale))
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


def main():
    args = parse_args()

    with args.coordinates.open() as coords_file:
        coordinates = parse_coordinates(coords_file)

    with args.operons.open() as operons_file:
        operons_df = parse_list_of_operons(operons_file, coordinates)

    gene_coordinates = find_oantigenes_data_from_gff(args.annotation)
    selected_operons = find_oantigene_operon_numbers(operons_df, gene_coordinates)
    selected_operons = sorted(selected_operons)

    records = gff_records(args.annotation)

    n = len(selected_operons)
    n2 = n // 2 + 1

    axes = [None] * n
    if args.plot:
        _, axes = plt.subplots(
            nrows=n2,
            ncols=2,
            figsize=(14, 2 * n2),
            layout="constrained",
        )
        axes = axes.flatten()

    for ax, operon_number in zip(axes, selected_operons):
        start, end, contig = get_operon_boundaries(operons_df, operon_number)

        contig_key = gff_key(contig)
        record = records[contig_key]

        # Do not crop the start of the gene
        start = start - 1
        fig = draw_region_by_coordinates(
            ax,
            operon_number,
            start,
            end,
            record,
            args,
        )

        fig.savefig(args.out_dir / f"operon_{operon_number}__{record.id}.png")
        plt.close(fig)

    if args.plot:
        plt.show()


if __name__ == "__main__":
    main()
