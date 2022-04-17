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
    container:
        containers["fastp"]
    params:
        opts="-l 30 --cut_tail --dedup --stdout",
        script=scripts_dir / "illumina_preprocess.sh",
    shadow:
        "shallow"
    shell:
        """
        bash {params.script} -r {wildcards.run} -i {input.run_info} -R {output.report} \
            -o {output.fastq} -t {threads} {params.opts} 2> {log}
        """


# TODO: use -p option in bwa mem for smart pairing https://manpages.org/bwa
