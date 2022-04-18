import sys
from typing import Dict, List

sys.stderr = open(snakemake.log[0], "w")

import json


def eprint(msg):
    print(msg, file=sys.stderr)


with open(snakemake.input.hunt2019) as fp:
    hunt_data = json.load(fp)

with open(snakemake.input.who2021) as fp:
    who_data = json.load(fp)

hunt_vars = set(hunt_data.keys())
eprint(f"{len(hunt_vars)} variants in Hunt2019 panel")
who_vars = set(who_data.keys())
eprint(f"{len(who_vars)} variants in WHO2021 panel")

common_vars = hunt_vars.intersection(who_vars)

eprint(f"{len(common_vars)} variants in common between panels")

var2drug: Dict[str, List[str]] = dict()

for v in common_vars:
    hunt_drugs = hunt_data[v]
    who_drugs = who_data[v]

    drugs = list(set(hunt_drugs + who_drugs))
    var2drug[v] = drugs

hunt_only_vars = hunt_vars - common_vars
eprint(f"{len(hunt_only_vars)} variants exclusive to Hunt2019 panel")
who_only_vars = who_vars - common_vars
eprint(f"{len(who_only_vars)} variants exclusive to WHO2021 panel")

for v in hunt_only_vars:
    var2drug[v] = hunt_data[v]

for v in who_only_vars:
    var2drug[v] = who_data[v]


eprint(f"{len(var2drug)} variants in combined panel")

with open(snakemake.output.resistance_json, "w") as fp:
    json.dump(var2drug, fp, indent=4, sort_keys=True)
