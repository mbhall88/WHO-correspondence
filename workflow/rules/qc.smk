rule illumina_preprocessing:
    input:
        run_dir=rules.download_data.output.outdir,
        run_info=rules.aggregate_run_info.output.run_info,
    output:
        fastq=results / "preprocessing/{run}/{run}.fq.gz",
        report=results / "preprocessing/{run}/{run}.report.html",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(4 * GB),
    log:
        log_dir / "illumina_preprocessing/{run}.log",
    conda:
        str(env_dir / "fastp.yaml")
    params:
        opts="-z 6 -l 30 --cut_tail --dedup --stdout",
        indelim=rules.aggregate_run_info.params.delim,
    shadow:
        "shallow"
    script:
        str(scripts_dir / "fastp.py")


# TODO: use -p option in bwa mem for smart pairing https://manpages.org/bwa
