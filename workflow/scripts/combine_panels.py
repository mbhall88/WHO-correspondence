import sys

sys.stderr = open(snakemake.log[0], "w")


def eprint(msg):
    print(msg, file=sys.stderr)


hunt_panel = set()

with open(snakemake.input.hunt2019) as fp:
    for line in map(str.rstrip, fp):
        gene, mut, alpha = line.split("\t")
        hunt_panel.add((gene, mut, alpha))

eprint(f"{len(hunt_panel)} mutations in the Hunt2019 panel")

who_panel = set()
with open(snakemake.input.who2021) as fp:
    for line in map(str.rstrip, fp):
        gene, mut, alpha = line.split("\t")
        who_panel.add((gene, mut, alpha))

eprint(f"{len(who_panel)} mutations in the WHO2021 panel")

common_vars = hunt_panel.intersection(who_panel)
eprint(f"{len(common_vars)} mutations in common between the panels")

combined_panel = who_panel.union(hunt_panel)

eprint(f"{len(combined_panel)} mutations in combined panel")

with open(snakemake.output.panel, "w") as fp:
    for t in combined_panel:
        print("\t".join(t), file=fp)
