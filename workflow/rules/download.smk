rule download_data:
    output:
        outdir=directory(data_dir / "fastq/{proj}/{sample}/{run}"),
        run_info=data_dir / "fastq/{proj}/{sample}/{run}/fastq-run-info.json",
    log:
        log_dir / "download_data/{proj}/{sample}/{run}.log",
    container:
        containers["fastq_dl"]
    resources:
        mem_mb=int(0.5 * GB),
    params:
        db="ena",
    shell:
        "fastq-dl --outdir {output.outdir} {wildcards.run} {params.db} > {log} 2>&1"


rule aggregate_run_info:
    input:
        dirs=expand(
            data_dir / "fastq/{proj}/{sample}/{run}",
            zip,
            run=RUNS,
            sample=SAMPLES,
            proj=PROJECTS,
        ),
    output:
        run_info=data_dir / "run_info.tsv",
    log:
        log_dir / "aggregate_run_info.log",
    params:
        delim="\t",
    script:
        str(scripts_dir / "aggregate_run_info.py")
