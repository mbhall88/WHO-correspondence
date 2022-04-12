rule illumina_preprocessing:
    input:
        run_dir=rules.download_data.output.outdir,
        run_info=rules.aggregate_run_info.output.run_info,
    output:
        fastq=directory(results / "preprocessing/{run}/{run}.fq.gz"),
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(4 * GB),
    log:
        log_dir / "illumina_preprocessing/{run}.log",
    container:
        containers["fastp"]
    params:
        opts="-z 6 -l 30 --cut_tail",
        indelim=rules.aggregate_run_info.params.delim
    shadow:
        "shallow"
    run:
        from pathlib import Path
        with open(input.run_info) as fp:
            found_this_run = False
            for line in map(str.rstrip, fp):
                run, layout = line.split(params.indelim)
                if run == wildcards.run:
                    found_this_run = True
                    break

        if not found_this_run:
            raise ValueError("Didn't find run in run_info")

        d = Path(input.run_dir)
        opts = params.opts
        inputs = "-i "
        outputs = ""
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
            opts += " --merge"
        else:
            raise KeyError(f"Got unknown library layout {layout}")

        # todo https://github.com/OpenGene/fastp#merge-paired-end-reads to use the correct output flags

        #
        # cmd="""
        # fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} \
        #   -w {threads} -h {output.report} 2> {log}
        # """

# todo: use -p option in bwa mem for smart pairing https://manpages.org/bwa