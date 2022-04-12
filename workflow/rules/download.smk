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


rule aggregate_run_info:
    input:
        dirs=expand(data_dir / "fastq/{run}", run=list(samplesheet.index))
    output:
        run_info=data_dir / "run_info.tsv"
    log:
        log_dir / "aggregate_run_info.log"
    params:
        delim="\t"
    script:
        str(scripts_dir / "aggregate_run_info.py")