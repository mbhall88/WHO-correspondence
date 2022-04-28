rule report:
    input:
        qc=rules.qc_summary.output.summary,
        csvs=expand(results / "amr_predictions/{panel}.csv", panel=PANELS),
        samplesheet=config["samplesheet"],
    output:
        table=results / "report/results.csv",
        plots=multiext(str(results / "report/plot"), ".png", ".svg"),
    resources:
        mem_mb=2 * GB,
    params:
        figsize=(13, 8),
        dpi=300,
        cov_threshold=config["cov_threshold"],
        ignore_drugs={"Ciprofloxacin", "Ofloxacin"},
        conf_interval=config["conf_interval"],
        minor_is_susceptible=config["minor_is_susceptible"],
        panel_names={
            "who2021": "WHO only",
            "hunt2019": "Mykrobe",
            "hall2022": "Combined",
        },
    log:
        notebook=results / "report/results.ipynb",
    conda:
        str(env_dir / "report.yaml")
    notebook:
        str(notebooks / "analysis.py.ipynb")
