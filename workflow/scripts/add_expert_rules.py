"""The WHO catalogue does not explicitly include 'expert rules'. This script adds those
rules to the catalogue.

The expert rules are:
- Any rpoB RRDR mutation (except synonymous mutations) for rifampicin. Here RRDR is
defined as rpoB codons 426-452
- Any premature stop codons or indels in the coding region of katG (isoniazid), ethA
(ethionamide), gid (streptomycin) and pncA (pyrazinamide)

From the paper, these expert rules were considered Group 2 (Associated with R – Interim)
"""
import sys
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from itertools import repeat
from pathlib import Path
from typing import Optional, Dict, TextIO, List, Tuple

import click
from loguru import logger


LOG_FMT = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)
Contig = str
Seq = str
Index = Dict[Contig, Seq]

codon2amino = {
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGC": "C",
    "TGT": "C",
    "TGA": "*",
    "TGG": "W",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
}


class Strand(Enum):
    Forward = "+"
    Reverse = "-"
    NotRelevant = "."
    Unknown = "?"

    def __str__(self) -> str:
        return str(self.value)


class RuleType(Enum):
    NonSynonymous = "nonsyn"
    Frameshift = "frame"
    Stop = "stop"


@dataclass
class Rule:
    rule_type: RuleType
    gene: str
    drugs: List[str]
    start: Optional[int] = None
    stop: Optional[int] = None
    grade: Optional[int] = None


def translate(seq: str, stop_last=True) -> str:
    if len(seq) % 3 != 0:
        raise ValueError("Sequence length must be a multiple of 3")

    prot = ""
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        prot += codon2amino[codon]

    if stop_last and not prot.endswith("*"):
        raise ValueError("Sequence did not end in a stop codon")

    return prot


@dataclass
class GffFeature:
    seqid: Contig
    source: str
    method: str  # correct term is type, but that is a python reserved variable name
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    score: float
    strand: Strand
    phase: int
    attributes: Dict[str, str]

    @staticmethod
    def from_str(s: str) -> "GffFeature":
        fields = s.split("\t")
        score = 0 if fields[5] == "." else float(fields[5])
        phase = -1 if fields[7] == "." else int(fields[7])
        attr_fields = fields[-1].split(";")
        attributes = {k: v for k, v in map(str.split, attr_fields, repeat("="))}
        return GffFeature(
            seqid=fields[0],
            source=fields[1],
            method=fields[2],
            start=int(fields[3]),
            end=int(fields[4]),
            score=score,
            strand=Strand(fields[6]),
            phase=phase,
            attributes=attributes,
        )

    def slice(self, zero_based: bool = True) -> Tuple[int, int]:
        """Get a tuple for slicing a python object.
        The reason this method is required is that GFF uses 1-based INCLUSIVE
        coordinates. Meaning the end position is also included in the slice.
        """
        if zero_based:
            return self.start - 1, self.end
        return self.start, self.end + 1

    def _extract_sequence(self, index: Index) -> str:
        pass

    def nucleotide_sequence(self, index: Index) -> str:
        return self._extract_sequence(index)

    def protein_sequence(self, index: Index) -> str:
        nuc_seq = self.nucleotide_sequence(index)
        return translate(nuc_seq)


def is_header(s: str) -> bool:
    if not s:
        return False
    return s[0] == ">"


class DuplicateContigsError(Exception):
    pass


def index_fasta(stream: TextIO) -> Index:
    fasta_index: Index = dict()
    sequence: List[Seq] = []
    name: Contig = ""
    for line in map(str.rstrip, stream):
        if not line:
            continue
        if is_header(line):
            if sequence and name:
                fasta_index[name] = "".join(sequence)
                sequence = []
            name = line.split()[0][1:]
            if name in fasta_index:
                raise DuplicateContigsError(
                    f"Contig {name} occurs multiple times in the fasta file."
                )
            continue
        else:
            sequence.append(line)
    if name and sequence:
        fasta_index[name] = "".join(sequence)

    return fasta_index


def setup_logging(verbose: bool) -> None:
    log_lvl = "INFO"
    if verbose:
        log_lvl = "DEBUG"
    logger.remove()
    logger.add(sys.stderr, level=log_lvl, format=LOG_FMT)


@click.command()
@click.option(
    "-a",
    "--annotation",
    help="GFF file with gene annotations",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "-r",
    "--reference",
    help="Reference FASTA",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--header/--no-header",
    default=True,
    help="Whether to include a header in the output",
)
@click.option(
    "-x",
    "--rules",
    help=(
        "Comma-separated file with expert rules. The format of this file is type, gene, start, "
        "stop, drugs (semi-colon (;) separated), grade. Valid types are nonsyn "
        "(non-synonymous mutations), frame (any frameshift "
        "indel), or stop (stop codon). If both start and stop are empty, the whole "
        "gene is used. If only start is given, then stop is considered the end of "
        "the gene and vice versa. Start and stop are CODONS, not positions, and are "
        "both inclusive. Grade is the (optional) grading to provide mutations arising from the rule"
    ),
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "-o",
    "--output",
    help="File to write output to",
    default="-",
    type=click.File(mode="w"),
)
@click.option("-v", "--verbose", help="Turns on debug-level logger.", is_flag=True)
@click.help_option("--help", "-h")
def main(
    rules: Path,
    annotation: Path,
    reference: Path,
    verbose: bool,
    header: bool,
    output: TextIO,
):
    setup_logging(verbose)

    logger.info("Indexing reference FASTA...")
    with reference.open() as ref_fp:
        index = index_fasta(ref_fp)
    logger.success(f"{len(index)} contig(s) indexed in the input file.")

    logger.info("Loading expert rules...")
    gene2rules = defaultdict(list)
    n_rules = 0
    with rules.open() as rules_fp:
        for row in map(str.rstrip, rules_fp):
            if row.startswith("type,"):  # header
                continue
            fields = row.split(",")
            start = int(fields[2]) if fields[2] else None
            stop = int(fields[3]) if fields[3] else None
            grade = int(fields[5]) if fields[5] else None
            drugs = ";".split(fields[4])
            rule = Rule(
                rule_type=RuleType(fields[0]),
                gene=fields[1],
                drugs=drugs,
                start=start,
                stop=stop,
                grade=grade,
            )
            gene2rules[rule.gene].append(rule)
            n_rules += 1

    logger.success(f"Loaded {n_rules} expert rules for {len(gene2rules)} genes")

    features = dict()
    with annotation.open() as gff_fp:
        for line in map(str.rstrip, gff_fp):
            if not line or line.startswith("#"):
                continue

            feature = GffFeature.from_str(line)
            if feature.method != "gene":
                continue

            if "Name" in feature.attributes:
                name = feature.attributes["Name"]
            elif "ID" in feature.attributes:
                name = feature.attributes["ID"]
            elif "gene" in feature.attributes:
                name = feature.attributes["gene"]
            else:
                continue

            if name not in gene2rules:
                continue

            features[name] = feature


if __name__ == "__main__":
    main()
