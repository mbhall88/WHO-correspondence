import sys

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path

with open(snakemake.input.run_info) as fp:
    found_this_run = False
    for line in map(str.rstrip, fp):
        run, layout = line.split(snakemake.params.indelim)
        if run == snakemake.wildcards.run:
            found_this_run = True
            break

if not found_this_run:
    raise ValueError("Didn't find run in run_info")

d = Path(snakemake.input.run_dir)
opts = snakemake.params.opts
inputs = "-i "
if layout == "SINGLE":
    fq = d / f"{run}.fastq.gz"
    if not fq.exists():
        raise FileNotFoundError(f"Expected single-end library {fq}")
    inputs += f"{fq}"
elif layout == "PAIRED":
    r1 = d / f"{run}_1.fastq.gz"
    if not r1.exists():
        raise FileNotFoundError(f"Expected R1 file {r1}")
    r2 = d / f"{run}_2.fastq.gz"
    if not r2.exists():
        raise FileNotFoundError(f"Expected R2 file {r2}")
    inputs += f"{r1} -I {r2}"
    opts += " --detect_adapter_for_pe"
else:
    raise KeyError(f"Got unknown library layout {layout}")

snakemake.shell(
    "( fastp -h {output.report} -w {threads} {inputs} {opts} | gzip -c ) > {output.fastq} 2> {log}"
)
