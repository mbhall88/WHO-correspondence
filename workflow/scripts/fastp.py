import sys

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path
import subprocess

with open(snakemake.input.run_info) as fp:
    found_this_run = False
    for line in map(str.rstrip, fp):
        run, files_str = line.split(snakemake.params.indelim)
        files = files_str.split(";")
        if run == snakemake.wildcards.run:
            found_this_run = True
            break

if not found_this_run:
    raise ValueError("Didn't find run in run_info")

d = Path(snakemake.input.run_dir)
opts = snakemake.params.opts
inputs = "-i "
if len(files) == 1:
    fq = Path(files[0])
    if not fq.exists():
        raise FileNotFoundError(f"Expected SE library {fq}")
    inputs += f"{fq}"
elif len(files) == 2:
    r1 = Path(files[0])
    if not r1.exists():
        raise FileNotFoundError(f"Expected R1 file {r1}")
    r2 = Path(files[1])
    if not r2.exists():
        raise FileNotFoundError(f"Expected R2 file {r2}")
    inputs += f"{r1} -I {r2}"
    opts += " --detect_adapter_for_pe"
else:
    raise KeyError(f"Got unknown library layout {files}")

fastp_cmd = f"fastp -h {snakemake.output.report} -w {snakemake.threads} {inputs} {opts}"

completed_proc = subprocess.run(
    f"( {fastp_cmd} | gzip -c ) > {snakemake.output.fastq}",
    shell=True,
    text=True,
    stderr=sys.stderr,
)

if completed_proc.returncode != 0:
    raise ChildProcessError("fastp command failed. Error should be above")
