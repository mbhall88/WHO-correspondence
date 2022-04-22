import sys

sys.stderr = open(snakemake.log[0], "w")

import json
import re
from collections import defaultdict


def split_var_name(name):
    items = re.match(r"([A-Z]+)([-0-9]+)([A-Z/\*]+)", name, re.I).groups()
    return items[0], int(items[1]), items[2]


with open(snakemake.input.resistance_json) as fp:
    var2drug = json.load(fp)


panel = defaultdict(list)
with open(snakemake.input.panel) as fp:
    for line in map(str.rstrip, fp):
        gene, mut, alpha = line.split("\t")
        ref, pos, alt = split_var_name(mut)
        panel[(gene, pos, alpha)].append(mut)

duplicates = set()
for (gene, pos, alpha), muts in panel.items():
    if len(muts) > 1 and any("X" in m for m in muts):
        for var in muts:
            ref, _, alt = split_var_name(var)
            if alt in ("X", "*"):
                continue
            duplicates.add((gene, var, alpha))
            full_var = f"{gene}_{var}"
            drugs = list(set(var2drug[full_var]))
            X_var = f"{gene}_{ref}{pos}X"
            var2drug[X_var].extend(drugs)
            var2drug[X_var] = list(sorted(set(var2drug[X_var])))
            del var2drug[full_var]


with open(snakemake.input.panel) as in_fp, open(snakemake.output.panel, "w") as out_fp:
    for line in map(str.rstrip, in_fp):
        gene, mut, alpha = line.split("\t")
        if (gene, mut, alpha) in duplicates:
            continue
        print(line, file=out_fp)

with open(snakemake.output.resistance_json, "w") as out_fp:
    json.dump(var2drug, out_fp, indent=4, sort_keys=True)
