rule download_data:
    output:
        outdir=directory(data_dir / "fastq/{run}"),
    log:
        log_dir / "download_data/{run}.log",
    container:
        containers["fastq_dl"]
    shadow:
        "shallow"
    resources:
        mem_mb=int(0.5 * GB),
    params:
        db="ena",
    shell:
        "fastq-dl --outdir {output.outdir} {wildcards.run} {params.db} > {log} 2>&1"
