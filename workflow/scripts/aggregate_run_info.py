import sys

sys.stderr = open(snakemake.log[0], "w")

import json
from pathlib import Path

DELIM = snakemake.params.delim


def eprint(msg):
    print(msg, file=sys.stderr)


out_data = []

for d in map(Path, snakemake.input.dirs):
    run = d.parts[-1]
    p = d / "fastq-run-info.json"
    with open(p) as fp:
        data = json.load(fp)[0]
        assert data["instrument_platform"] == "ILLUMINA", p

        layout = data["library_layout"]

        if layout == "SINGLE":
            fq = d / f"{run}.fastq.gz"
            if not fq.exists():
                raise FileNotFoundError(f"Expected single-end library {fq}")
        elif layout == "PAIRED":
            r1 = d / f"{run}_1.fastq.gz"
            if not r1.exists():
                raise FileNotFoundError(f"Expected R1 file {r1}")
            r2 = d / f"{run}_2.fastq.gz"
            if not r2.exists():
                raise FileNotFoundError(f"Expected R2 file {r2}")
            fq = d / f"{run}.fastq.gz"
            if fq.exists():
                eprint(f"[INFO]: Removing orphan SE file for paired-end accession {fq}")
        else:
            raise KeyError(f"Got unknown library layout {layout}")

        if data["tax_id"] != "1773":
            eprint(f"[WARNING]: Got non-MTB tax ID for {run} - {data['tax_id']}")

    out_data.append((run, layout))

with open(snakemake.output.run_info, "w") as fp:
    print(f"run{DELIM}layout", file=fp)
    for t in out_data:
        print(DELIM.join(t), file=fp)
