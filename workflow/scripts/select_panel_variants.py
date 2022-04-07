import sys

sys.stderr = open(snakemake.log[0], "w")
from collections import defaultdict
import json


def eprint(msg):
    print(msg, file=sys.stderr)


max_grade = int(snakemake.wildcards.grade)

var2drug = defaultdict(list)

n_kept = 0
n_total = 0
eprint(
    f"Filtering full panel and only keeping variants with a grading <= {max_grade}..."
)

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

        variant = f"{gene}_{mut}"
        var2drug[variant].append(drug.capitalize())

        print("\t".join([gene, mut, alpha]), file=out_fp)
        n_kept += 1

eprint(f"Kept {n_kept} of {n_total} variants")

eprint("Write variant to drug JSON file...")
with open(snakemake.output.resistance_json, "w") as fp:
    json.dump(var2drug, fp, indent=4, sort_keys=True)

eprint("Done!")
