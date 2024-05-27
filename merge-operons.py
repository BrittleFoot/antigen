import argparse
import io
import re
import sys
from dataclasses import dataclass
from functools import wraps
from pathlib import Path
from typing import Dict, Literal, NamedTuple

from BCBio import GFF

O_ANTIGEN_GENES = {
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
    ######
    "yibK",
    "wzc",
    "wzb",
    "wza",
    "galE",
    "gne",
    "wpaD",
    "ugd",
    "wpaC",
    "wpaB",
    "wzy",
    "wpaA",
    "wzx",
    "qdtB",
    "qdtA",
    "rmlA",
    "qdtf",
    "cpxA",
    "wec",
    "rffG",
    "rffH",
}

H_ANTIGEN_GENES = {
    "fliC",
    "flkA",
    "fllA",
    "flmA",
    "flgA",
    "flgB",
    "flgC",
    "flgD",
    "flgE",
    "flgF",
    "flhA",
    "flhB",
    "fliA",
    "fliB",
    "fliD",
    "fliE",
}

K_ANTIGEN_GENES = {
    "kpsM",
    "kpsT",
    "kpsD",
    "kpsE",
    "cpsA",
    "cpsB",
    "cpsC",
    "cpsD",
}

ANTIGEN_GENES = {
    "O": O_ANTIGEN_GENES,
    "H": H_ANTIGEN_GENES,
    "K": K_ANTIGEN_GENES,
}


@wraps(print)
def err(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


AntigenType = Literal["O", "H", "K"]


@dataclass
class InputArgs:
    type: AntigenType
    #
    annotation: Path
    operons: Path
    coordinates: Path
    #
    output: Path


@dataclass(frozen=True, slots=True)
class OperonGene:
    operon: int
    id_gene: str


NODE_ID_PATTERN = re.compile(r"(\d+)")


def print_feature(feature):
    q = feature.qualifiers
    loc = feature.location
    print(
        loc.start,
        loc.end,
        loc.strand,
        q.get("gene")[0],
        q.get("Name")[0],
        sep="\t",
    )


def contig_key(record):
    """parse contig id

    ie when we have NODE_36_length_1620_cov_15 from operonfinder, but
    contig_36 from annotation, parse 36 as key
    """
    ids = NODE_ID_PATTERN.findall(record.id)
    if ids:
        return int(ids[0])

    raise ValueError(f"No obvious id for {record}")


def mark_operonfinder_operons(records, operons: Dict[str, int]):
    for record in records:
        for feature in record.features:
            id = feature.qualifiers["ID"][0]
            if id not in operons:
                continue
            feature.qualifiers["operon"] = [operons[id]]

        yield record


def parse_list_of_operons(list_file: io.FileIO):
    skipheader = False
    current_operon = None

    for line in list_file:
        line = line.strip()

        if not skipheader:
            skipheader = True
            continue

        line = line.split("\t")

        if len(line) == 1:
            current_operon = int(line[0])
            continue

        id_gene = line[0]
        if id_gene == "NA":
            continue

        assert len(line) == 7, line

        yield OperonGene(current_operon, id_gene)


def is_antigen(feature, antigen_type: AntigenType):
    q = feature.qualifiers
    gene_names = (q.get("gene", []) + q.get("Gene", [])) or q.get("Name", [])

    if not gene_names:
        return False

    if len(gene_names) > 1:
        raise NotImplementedError(f"More than one gene name: {gene_names} in {feature}")

    gene = gene_names[0]

    # replace to just flat `in` if substring is not needed,
    # and the whole gene in 'gene' qualifier
    for antigene_name in ANTIGEN_GENES[antigen_type]:
        if antigene_name.find(gene) != -1:
            return True

    return False


def mark_antigenes(records, antigen_type: AntigenType):
    for record in records:
        for feature in record.features:
            antigen = is_antigen(feature, antigen_type=antigen_type)
            feature.qualifiers["antigen"] = [antigen]
            if antigen:
                feature.qualifiers["antigen_type"] = [antigen_type]

            if antigen:
                print(record.id, end="\t")
                print_feature(feature)

        yield record


def parse_gff(gff_file: io.FileIO):
    return GFF.parse(gff_file)


class Location(NamedTuple):
    start: int
    end: int
    strand: str

    @classmethod
    def from_feature(cls, feature):
        return cls(
            start=feature.location.start,
            end=feature.location.end,
            strand=feature.location.strand,
        )


def pairwise_gffs(a_records, b_records):
    b_keys = {contig_key(b): b for b in b_records}
    for a in a_records:
        yield a, b_keys.get(contig_key(a), None)


def match_gffs(annoation, operons):
    for ann, oper in pairwise_gffs(annoation, operons):
        if oper is not None:
            match_records(ann, oper)
        yield ann


def match_records(annotation, ofinder):
    loc_operons = {
        Location.from_feature(f): f.qualifiers["operon"][0]
        for f in ofinder.features
        if "operon" in f.qualifiers
    }

    for feature in annotation.features:
        loc = Location.from_feature(feature)

        if loc in loc_operons:
            feature.qualifiers["operon"] = [loc_operons[loc]]
            continue


def main(args: InputArgs):
    print(f"Looking for {args.type}-antigens in the annotation")
    with args.operons.open() as operonsh:
        operon_ids = {o.id_gene: o.operon for o in parse_list_of_operons(operonsh)}

    with args.annotation.open() as inputh, args.coordinates.open() as coordsh:
        gff = mark_antigenes(parse_gff(inputh), args.type)
        operons = mark_operonfinder_operons(parse_gff(coordsh), operon_ids)

        matched = match_gffs(gff, operons)

        with args.output.open("w") as outh:
            GFF.write(matched, out_handle=outh)


def parse_args(args: list[str] = None):
    arp = argparse.ArgumentParser(
        "merge-operons",
        description=(
            "combines annotation and the output of the operon finder into one gff"
            + "where each feature have an operon number and antigen mark qualifiers"
        ),
    )

    arp.add_argument(
        "-t",
        "--type",
        choices=["O", "H", "K"],
        required=True,
        help="Antigen type",
    )

    g = arp.add_argument_group("annotation file")
    g.add_argument("annotation", type=Path, help="gff annotation file path")

    g = arp.add_argument_group(
        "Operon mapper",
        description="The files you will get as a result from the Operon Mapper run",
    )
    g.add_argument("operons", type=Path, help="list_of_operons file path")
    g.add_argument("coordinates", type=Path, help="ORFs_coordinates file path")

    g = arp.add_argument_group("output options")
    g.add_argument("-o", "--output", type=Path, help="output file path")

    ###########################################################

    parsed = InputArgs(**arp.parse_args(args=args).__dict__)
    if parsed.output is None:
        prefix = "operons-"

        out_name = prefix + parsed.annotation.name
        parsed.output = parsed.annotation.parent / out_name

    return parsed


if __name__ == "__main__":
    main(parse_args())
