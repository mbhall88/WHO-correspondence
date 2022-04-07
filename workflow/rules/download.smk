rule download_data:
    output:
        r1=data_dir / "fastq/{run}/{run}_1.fastq.gz",
        r2=data_dir / "fastq/{run}/{run}_2.fastq.gz",
        run_info=data_dir / "fastq/{run}/fastq-run-info.json",
    log:
        log_dir / "download_data/{run}/{run}.log",
    container:
        containers["fastq_dl"]
    resources:
        mem_mb=int(0.5 * GB),
    params:
        db="sra",
        outdir=lambda wildcards, output: Path(output.run_info).parent,
    shell:
        "fastq-dl --outdir {params.outdir} {wildcards.run} {params.db} > {log} 2>&1"
