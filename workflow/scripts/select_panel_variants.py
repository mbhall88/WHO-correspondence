import sys
from typing import Dict, Tuple, Set

sys.stderr = open(snakemake.log[0], "w")
from collections import defaultdict
import json
import re
from pathlib import Path


def split_var_name(name: str) -> Tuple[str, int, str]:
    items = re.match(r"([A-Z]+)([-0-9]+)([A-Z/\*]+)", name, re.I).groups()
    return items[0], int(items[1]), items[2]


TRANSLATE = str.maketrans("ATGC", "TACG")


def revcomp(s: str) -> str:
    return complement(s)[::-1]


def complement(s: str) -> str:
    return s.upper().translate(TRANSLATE)


def eprint(msg):
    print(msg, file=sys.stderr)


def load_strands(path: Path) -> Dict[str, str]:
    strands = dict()
    with open(path) as fp:
        for line in map(str.rstrip, fp):
            if line.startswith("#") or not line:
                continue
            fields = line.split("\t")
            ftr = fields[2]
            if ftr != "gene":
                continue
            extras = fields[8].split(";")
            gene = ""
            for e in extras:
                if e.startswith("Name=") or e.startswith("gene="):
                    gene = e.split("=")[-1]
                    break
            if not gene:
                raise ValueError(f"Couldn't find the gene name for row {line}")
            strands[gene] = fields[6]

    return strands


max_grade = snakemake.params.max_grade
gene_strands = load_strands(snakemake.input.features)
adjustments = snakemake.params.adjustments
var2drug = defaultdict(list)

n_kept = 0
n_total = 0
eprint(
    f"Filtering full panel and only keeping variants with a grading <= {max_grade}..."
)

# a counter to make sure we see all expected adjustments
seen_adjustments: Set[str] = set()

with open(snakemake.input.full_panel) as in_fp, open(
    snakemake.output.panel, "w"
) as out_fp:
    _ = next(in_fp)  # skip header
    for line in map(str.rstrip, in_fp):
        n_total += 1
        gene, mut, alpha, drug, grading = line.split("\t")
        grading = int(grading)

        if grading > max_grade:
            continue

        if f"{gene}_{mut}" in adjustments:
            seen_adjustments.add(f"{gene}_{mut}")
            variant = adjustments[f"{gene}_{mut}"]
            gene, mut = variant.split("_")
        else:
            # The position of DNA mutations (inside genes) on the rev strand point to the *end* of the
            # reference allele. See the data cleaning notebook for an explanation
            is_promotor_mut = "-" in mut
            if gene_strands[gene] == "-" and alpha == "DNA" and not is_promotor_mut:
                ref, pos, alt = split_var_name(mut)
                ref = revcomp(ref)
                alt = revcomp(alt)
                pos -= len(ref) - 1
                mut = f"{ref}{pos}{alt}"

            variant = f"{gene}_{mut}"

        var2drug[variant].append(drug.capitalize())

        print("\t".join([gene, mut, alpha]), file=out_fp)
        n_kept += 1

eprint(f"Kept {n_kept} of {n_total} variants")

expected_adjustments = set(adjustments.keys())
diff = expected_adjustments - seen_adjustments
if diff:
    raise KeyError(f"Didn't see the following expected adjustments:\n{diff}")

eprint("Write variant to drug JSON file...")
with open(snakemake.output.resistance_json, "w") as fp:
    json.dump(var2drug, fp, indent=4, sort_keys=True)

eprint("Done!")
