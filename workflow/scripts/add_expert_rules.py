"""The WHO catalogue does not explicitly include 'expert rules'. This script adds those
rules to the catalogue.

The expert rules are:
- Any rpoB RRDR mutation (except synonymous mutations) for rifampicin. Here RRDR is
defined as rpoB codons 426-452
- Any premature stop codons or indels in the coding region of katG (isoniazid), ethA
(ethionamide), gid (streptomycin) and pncA (pyrazinamide)

From the paper, these expert rules were considered Group 2 (Associated with R – Interim)
"""
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from itertools import repeat
from pathlib import Path
from typing import Optional, Dict, TextIO, List, Tuple, Set

import click
from loguru import logger

LOG_FMT = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)
TRANSLATE = str.maketrans("ATGC", "TACG")
STOP = "*"
MISSENSE = "X"
PROT = "PROT"
DNA = "DNA"
NUCLEOTIDES = ["A", "C", "G", "T"]
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
    "TAA": STOP,
    "TAG": STOP,
    "TGC": "C",
    "TGT": "C",
    "TGA": STOP,
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


def revcomp(s: str) -> str:
    return complement(s)[::-1]


def complement(s: str) -> str:
    return s.upper().translate(TRANSLATE)


class Strand(Enum):
    Forward = "+"
    Reverse = "-"
    NotRelevant = "."
    Unknown = "?"

    def __str__(self) -> str:
        return str(self.value)


class RuleType(Enum):
    Missense = "missense"
    Frameshift = "frame"
    Nonsense = "nonsense"


@dataclass
class Rule:
    rule_type: RuleType
    gene: str
    drugs: List[str]
    start: Optional[int] = None  # 1-based inclusive
    stop: Optional[int] = None  # 1-based inclusive
    grade: Optional[int] = None

    def _apply_nonsense(
        self, nucleotide_seq: str
    ) -> Set[Tuple[str, str, str, str, int]]:
        protein = translate(nucleotide_seq, stop_last=True)
        start = 0 if self.start is None else self.start - 1
        stop = len(protein) if self.stop is None else self.stop
        grade = -1 if self.grade is None else self.grade
        mutations = set()
        for i in range(start, stop):
            aa = protein[i]
            if aa == STOP:
                continue
            mut = f"{aa}{i+1}{STOP}"
            for drug in self.drugs:
                mutations.add((self.gene, mut, PROT, drug, grade))

        return mutations

    def _apply_frameshift(
        self, nucleotide_seq: str
    ) -> Set[Tuple[str, str, str, str, int]]:
        """There is an assumption in this method that 1bp up and downstream of the gene
        have been included in nucleotide_seq
        """
        base_before = nucleotide_seq[0]
        base_after = nucleotide_seq[-1]
        nuc_seq = nucleotide_seq[1:-1]
        start_at = 0 if self.start is None else self.start * 3 - 3
        end_at = len(nuc_seq) if self.stop is None else self.stop * 3 - 3

        indels = get_indel_combinations(
            nuc_seq,
            base_before=base_before,
            base_after=base_after,
            start_at=start_at,
            end_at=end_at,
        )
        mutations = set()
        grade = -1 if self.grade is None else self.grade
        for drug in self.drugs:
            for ref, pos, alt in indels:
                mut = f"{ref}{pos}{alt}"
                mutations.add((self.gene, mut, DNA, drug, grade))

        return mutations

    def _apply_missense(
        self, nucleotide_seq: str
    ) -> Set[Tuple[str, str, str, str, int]]:
        return {
            (gene, mut.replace(STOP, MISSENSE), alpha, drug, grade)
            for gene, mut, alpha, drug, grade in self._apply_nonsense(nucleotide_seq)
        }

    def apply(self, nucleotide_seq: str) -> Set[Tuple[str, str, str, str, int]]:
        if self.rule_type is RuleType.Nonsense:
            return self._apply_nonsense(nucleotide_seq)
        elif self.rule_type is RuleType.Frameshift:
            return self._apply_frameshift(nucleotide_seq)
        elif self.rule_type is RuleType.Missense:
            return self._apply_missense(nucleotide_seq)
        else:
            raise NotImplementedError(f"Don't know how to apply {self.rule_type}")


def get_indel_combinations(
    seq: str,
    base_before: str,
    base_after: str,
    start_at: int = 0,
    end_at: Optional[int] = None,
) -> List[Tuple[str, int, str]]:
    indels = []
    if end_at is None:
        end_at = len(seq)

    for i in range(start_at, end_at):
        pos = i
        for indel_length in [1, 2]:
            if i == 0:
                ref = base_before + seq[i : i + indel_length]
                pos = -1
            elif i == len(seq) - 1 and indel_length == 2:
                ref = seq[i - 1 : i + indel_length] + base_after
            else:
                ref = seq[i - 1 : i + indel_length]
            alt = ref[0]
            indels.append((ref, pos, alt))
            if indel_length == 1:
                continue

            for nuc1 in NUCLEOTIDES:
                indels.append((alt, pos, alt + nuc1))
                for nuc2 in NUCLEOTIDES:
                    indels.append((alt, pos, alt + nuc1 + nuc2))
    return indels


def split_var_name(name: str) -> Tuple[str, int, str]:
    items = re.match(r"([A-Z]+)([-0-9]+)([A-Z/\*]+)", name, re.I).groups()
    return items[0], int(items[1]), items[2]


def translate(seq: str, stop_last=True) -> str:
    if len(seq) % 3 != 0:
        raise ValueError("Sequence length must be a multiple of 3")

    prot = ""
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        prot += codon2amino[codon]

    if stop_last and not prot.endswith(STOP):
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

    def _extract_sequence(
        self, index: Index, start_offset: int = 0, end_offset: int = 0
    ) -> str:
        refseq = index.get(self.seqid)
        if refseq is None:
            raise IndexError(f"Contig {self.seqid} does not exist in reference")
        s, e = self.slice(zero_based=True)
        s -= start_offset
        e += end_offset
        return refseq[s:e]

    def nucleotide_sequence(
        self, index: Index, start_offset: int = 0, end_offset: int = 0
    ) -> str:
        nuc_seq = self._extract_sequence(
            index, start_offset=start_offset, end_offset=end_offset
        )
        if self.strand is Strand.Reverse:
            nuc_seq = revcomp(nuc_seq)

        return nuc_seq

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
        "stop, drugs (semi-colon (;) separated), grade. Valid types are missense, "
        "frame (any frameshift indel), or nonsense (premature stop codon). If both "
        "start and stop are empty, the whole "
        "gene is used. If only start is given, then stop is considered the end of "
        "the gene and vice versa. Start and stop are CODONS, not positions, and are "
        "both 1-based inclusive. Grade is the (optional) grading to provide mutations arising from the rule"
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
@click.option("-d", "--delim", help="Output file delimiter", default="\t")
@click.option("-v", "--verbose", help="Turns on debug-level logger.", is_flag=True)
@click.help_option("--help", "-h")
def main(
    rules: Path,
    annotation: Path,
    reference: Path,
    verbose: bool,
    header: bool,
    delim: str,
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
            drugs = fields[4].split(";")
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

    panel = set()
    for gene, rules in gene2rules.items():
        ftr = features[gene]

        for rule in rules:
            if rule.rule_type is RuleType.Frameshift:
                nuc_seq = ftr.nucleotide_sequence(
                    index=index, start_offset=1, end_offset=1
                )
            else:
                nuc_seq = ftr.nucleotide_sequence(index)
            rule_variants = rule.apply(nuc_seq)
            if rule_variants:
                logger.debug(
                    f"Generated {len(rule_variants):,} variants for rule {rule}"
                )
                panel = panel.union(rule_variants)

    if header:
        print(
            delim.join(["gene", "mutation", "alphabet", "drug", "grading"]), file=output
        )

    for variant in sorted(panel, key=lambda t: (t[0], split_var_name(t[1])[1], t[1])):
        print(delim.join(map(str, variant)), file=output)

    logger.success(f"Generated {len(panel):,} variants in total")


if __name__ == "__main__":
    main()
