rule illumina_preprocessing:
    input:
        run_dir=rules.download_data.output.outdir,
        run_info=rules.aggregate_run_info.output.run_info,
    output:
        fastq=results / "preprocessing/{run}/{run}.fq.gz",
        report=results / "preprocessing/{run}/{run}.report.html",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(8 * GB),
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


rule build_decontamination_db:
    output:
        fasta=decontam_db / "remove_contam.fa.gz",
        metadata=decontam_db / "remove_contam.tsv",
    threads: 1
    resources:
        mem_mb=GB,
    params:
        script=scripts_dir / "download_tb_reference_files.pl",
        outdir=lambda wildcards, output: Path(output.fasta).parent,
    conda:
        str(env_dir / "decontam_db.yaml")
    log:
        log_dir / "build_decontamination_db.log",
    shadow:
        "shallow"
    shell:
        "perl {params.script} {params.outdir} &> {log}"


BWA_EXTNS = [".amb", ".ann", ".bwt", ".pac", ".sa"]


rule index_decontam_db:
    input:
        fasta=rules.build_decontamination_db.output.fasta,
    output:
        index=multiext(str(decontam_db / "remove_contam.fa.gz"), *BWA_EXTNS),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(32 * GB),
    log:
        log_dir / "index_decontam_db.log",
    conda:
        str(env_dir / "aln_tools.yaml")
    shell:
        "bwa index {input.fasta} 2> {log}"


rule map_to_decontam_db:
    input:
        index=rules.index_decontam_db.output.index,
        ref=rules.index_decontam_db.input.fasta,
        reads=rules.illumina_preprocessing.output.fastq,
        run_info=rules.aggregate_run_info.output.run_info,
    output:
        bam=temp(results / "mapped/{run}/{run}.sorted.bam"),
        index=temp(results / "mapped/{run}/{run}.sorted.bam.bai"),
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(12 * GB),
    params:
        map_extras="-M",
        script=scripts_dir / "illumina_preprocess.sh",
    conda:
        str(env_dir / "aln_tools.yaml")
    log:
        log_dir / "map_to_decontam_db/{run}/{run}.log",
    shell:
        """
        bash {params.script} -r {wildcards.run} -i {input.run_info} -R {input.reads} \
            -o {output.bam} -d {input.ref} -t {threads} {params.map_extras} 2> {log}
        """


rule filter_contamination:
    input:
        bam=rules.map_to_decontam_db.output.bam,
        metadata=rules.build_decontamination_db.output.metadata,
    output:
        keep_ids=results / "filtered/{run}/keep.reads",
        contam_ids=results / "filtered/{run}/contaminant.reads",
        unmapped_ids=results / "filtered/{run}/unmapped.reads",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * GB,
    conda:
        str(env_dir / "filter.yaml")
    params:
        script=scripts_dir / "filter_contamination.py",
        extra="--ignore-secondary",
        outdir=lambda wildcards, output: Path(output.keep_ids).parent,
    log:
        log_dir / "filter_contamination/{run}/{run}.log",
    shell:
        """
        python {params.script} {params.extra} \
            -i {input.bam} \
            -m {input.metadata} \
            -o {params.outdir} 2> {log}
        """


rule extract_decontaminated_reads:
    input:
        reads=rules.map_to_decontam_db.input.reads,
        read_ids=rules.filter_contamination.output.keep_ids,
    output:
        reads=results / "filtered/{run}/{run}.filtered.fq.gz",
        stats=results / "filtered/{run}/{run}.filtered.stats.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: int(4 * GB) * attempt,
    log:
        log_dir / "extract_decontaminated_reads/{run}/{run}.log",
    container:
        containers["seqkit"]
    shell:
        """
        seqkit grep -o {output.reads} -f {input.read_ids} {input.reads} 2> {log}
        seqkit stats -a -T {output.reads} > {output.stats} 2>> {log}
        """
