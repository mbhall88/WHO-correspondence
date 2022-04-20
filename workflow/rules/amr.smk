rule mykrobe_who_panel:
    input:
        reads=results / "filtered/{proj}/{sample}/{run}/{run}.filtered.fq.gz",
        probes=panel_dir / "who2021/grading_2/probes.fa",
        resistance_json=panel_dir / "who2021/grading_2/var2drug.json",
    output:
        report=results / "amr_predictions/who2021/{proj}/{sample}/{run}.mykrobe.json",
    shadow:
        "shallow"
    resources:
        mem_mb=int(2 * GB),
    container:
        containers["mykrobe"]
    group:
        "mykrobe"
    log:
        log_dir / "mykrobe_who_panel/{proj}/{sample}/{run}.log",
    params:
        opts=" ".join(
            [
                "--force",
                "-O json",
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
        reads=results / "filtered/{proj}/{sample}/{run}/{run}.filtered.fq.gz",
    output:
        report=results / "amr_predictions/hunt2019/{proj}/{sample}/{run}.mykrobe.json",
    shadow:
        "shallow"
    resources:
        mem_mb=int(2 * GB),
    container:
        containers["mykrobe"]
    group:
        "mykrobe"
    log:
        log_dir / "mykrobe_default_panel/{proj}/{sample}/{run}.log",
    params:
        opts=" ".join(
            [
                "--force",
                "-O json",
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
        reads=results / "filtered/{proj}/{sample}/{run}/{run}.filtered.fq.gz",
        probes=rules.construct_combined_panel.output.probes,
        resistance_json=rules.construct_combined_panel.input.resistance_json,
    output:
        report=results / "amr_predictions/hall2022/{proj}/{sample}/{run}.mykrobe.json",
    shadow:
        "shallow"
    group:
        "mykrobe"
    resources:
        mem_mb=int(2 * GB),
    container:
        containers["mykrobe"]
    log:
        log_dir / "mykrobe_combined_panel/{proj}/{sample}/{run}.log",
    params:
        opts=" ".join(
            [
                "--force",
                "-O json",
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


rule combine_amr_reports:
    input:
        reports=expand(
            results / "amr_predictions/{{panel}}/{proj}/{sample}/{run}.mykrobe.json",
            run=RUNS,
            sample=SAMPLES,
            proj=PROJECTS,
        ),
    output:
        report=results / "amr_predictions/{panel}.csv",
    log:
        log_dir / "combine_amr_reports/{panel}.log"
    container:
        containers["python"]
    script:
        str(scripts_dir / "combine_amr_reports.py")
