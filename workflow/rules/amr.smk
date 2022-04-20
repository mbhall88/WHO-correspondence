rule mykrobe_who_panel:
    input:
        reads=rules.extract_decontaminated_reads.output.reads,
        probes=panel_dir / "who2021/grading_2/probes.fa",
        resistance_json=panel_dir / "who2021/grading_2/var2drug.json",
    output:
        report=results / "amr_predictions/who2021/{proj}/{sample}/{run}.mykrobe.tsv",
    shadow:
        "shallow"
    resources:
        mem_mb=int(2 * GB),
    container:
        containers["mykrobe"]
    group:
        "mykrobe"
    log:
        rule_log_dir / "mykrobe_who_panel/{proj}/{sample}/{run}.log",
    params:
        opts=" ".join(
            [
                "--force",
                "-O tsv",
                "-D 0.20",
                "--species custom",
                "-e 0.001",
                "--panel custom",
            ]
        ),
    threads: 2
    shell:
        """
        mykrobe predict {params.opts} -o {output.report} -i {input.reads} \
          --sample {wildcards.run} -t {threads} -m {resources.mem_mb}MB \
          -P {input.probes} -R {input.resistance_json} > {log} 2>&1
        """


rule mykrobe_default_panel:
    input:
        reads=rules.extract_decontaminated_reads.output.reads,
    output:
        report=results / "amr_predictions/hunt2019/{proj}/{sample}/{run}.mykrobe.tsv",
    shadow:
        "shallow"
    resources:
        mem_mb=int(2 * GB),
    container:
        containers["mykrobe"]
    group:
        "mykrobe"
    log:
        rule_log_dir / "mykrobe_default_panel/{proj}/{sample}/{run}.log",
    params:
        opts=" ".join(
            [
                "--force",
                "-O tsv",
                "-D 0.20",
                "--species tb",
                "-e 0.001",
            ]
        ),
    threads: 2
    shell:
        """
        mykrobe predict {params.opts} -o {output.report} -i {input.reads} \
          --sample {wildcards.run} -t {threads} -m {resources.mem_mb}MB > {log} 2>&1
        """


rule mykrobe_combined_panel:
    input:
        reads=rules.extract_decontaminated_reads.output.reads,
        probes=rules.construct_combined_panel.output.probes,
        resistance_json=rules.construct_combined_panel.input.resistance_json,
    output:
        report=results / "amr_predictions/hall2022/{proj}/{sample}/{run}.mykrobe.tsv",
    shadow:
        "shallow"
    group:
        "mykrobe"
    resources:
        mem_mb=int(2 * GB),
    container:
        containers["mykrobe"]
    log:
        rule_log_dir / "mykrobe_combined_panel/{proj}/{sample}/{run}.log",
    params:
        opts=" ".join(
            [
                "--force",
                "-O tsv",
                "-D 0.20",
                "--species custom",
                "-e 0.001",
                "--panel custom",
            ]
        ),
    threads: 2
    shell:
        """
        mykrobe predict {params.opts} -o {output.report} -i {input.reads} \
          --sample {wildcards.run} -t {threads} -m {resources.mem_mb}MB \
          -P {input.probes} -R {input.resistance_json} > {log} 2>&1
        """
