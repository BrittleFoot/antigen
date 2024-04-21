import argparse
import io
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict

from BCBio import GFF
from dna_features_viewer import BiopythonTranslator
from matplotlib import pyplot as plt


@dataclass
class InputArgs:
    annotation: Path
    #
    output: Path
    label: str
    no_colors: bool


def parse_gff(gff_file: io.FileIO):
    return GFF.parse(gff_file)


def qual(feature, qualifier_key, default=None):
    return feature.qualifiers.get(qualifier_key, [default])[0]


def draw_region(
    record,
    *,
    start: int,
    end: int,
    title: str,
    args: InputArgs,
):
    """
    Function produces schematic plot of DNA region specified with coordinates;
    visualizes only gene feature types and labels them
    """
    labels = [
        "label",
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
    translator.ignored_features_types = ["region", "remark"]
    translator.label_fields = labels
    graphic_record = translator.translate_record(record)

    operon = graphic_record.crop((start, end))
    operon.plot(elevate_outline_annotations=False)
    ax = plt.gca()
    ax.set_title(title)

    plt.savefig(args.output / f"{title}.png")


def group_by_operons(features):
    groups: Dict[str, Any] = defaultdict(list)

    for feature in features:
        if operon := qual(feature, "operon"):
            groups[operon].append(feature)

    return groups


def filter_by_oantigens(features):
    for feature in features:
        if qual(feature, "oantigen") == "True":
            return True

    return False


def oantigen_operons(operons):
    for operon, features in operons.items():
        if filter_by_oantigens(features):
            yield operon, features


def boundaries(features):
    start = None
    end = None

    for feature in features:
        loc = feature.location

        if start is None and end is None:
            start = loc.start
            end = loc.end

        start = min(start, loc.start)
        end = max(end, loc.end)

    return start, end


def all_oantigen_operons(records):
    for record in records:
        operon_groups = group_by_operons(record.features)
        for operon, features in oantigen_operons(operon_groups):
            yield record, operon, boundaries(features)


def colorize(records):
    for record in records:
        for feature in record.features:
            if qual(feature, "oantigen") == "True":
                feature.qualifiers["color"] = "#dc6678"
        yield record


def main(args: InputArgs):
    with args.annotation.open() as annotation:
        ann = parse_gff(annotation)

        if not args.no_colors:
            ann = colorize(ann)

        for record, operon, (start, end) in all_oantigen_operons(ann):
            draw_region(
                record,
                start=start,
                end=end,
                title=f"Contig: {record.id}, Operon: {operon}",
                args=args,
            )


def parse_args(args: list[str] = None):
    arp = argparse.ArgumentParser(
        "draw-operons",
        usage="Draws ",
    )

    g = arp.add_argument_group("Gff input")
    g.add_argument(
        "annotation",
        type=Path,
        help="gff annotation from the output of the merge-operons.py",
    )

    g = arp.add_argument_group("Output options")
    g.add_argument("-o", "--output", type=Path, help="output dir path")
    g.add_argument(
        "-l",
        "--label",
        default="gene",
        help="which label for the gene to use (from gff qualifiers)",
    )
    g.add_argument(
        "--no-colors",
        action="store_true",
        help="draw oantigen genes with outstanding colors",
    )

    ###########################################################

    parsed = InputArgs(**arp.parse_args(args=args).__dict__)
    if parsed.output is None:
        prefix = "pictures-"
        out_name = prefix + parsed.annotation.name
        parsed.output = parsed.annotation.parent / out_name
    parsed.output.mkdir(exist_ok=True, parents=True)

    return parsed


if __name__ == "__main__":
    main(parse_args())
