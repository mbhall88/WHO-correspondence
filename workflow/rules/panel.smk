"""Rules relating to the construction of the necessary Mykrobe panels"""


rule extract_background_vcfs:
    input:
        tarball=config["background_variants"],
    output:
        vcf_dir=directory(panel_dir / "background_vcfs"),
    log:
        log_dir / "extract_background_vcfs.log",
    container:
        containers["base"]
    shell:
        """
        mkdir -p {output} 2> {log}
        tar xzf {input} -C {output.vcf_dir} --strip-components 1 2>> {log}
        """


rule select_panel_variants:
    input:
        full_panel=full_panel,
        features=features
    output:
        panel=panel_dir / "who2021/grading_{grade}/panel.tsv",
        resistance_json=panel_dir / "who2021/grading_{grade}/var2drug.json",
    log:
        log_dir / "select_panel_variants/grade_{grade}.log",
    params:
        adjustments={
            "eis_CG-7C": "eis_CG-8C",
            "pncA_TC-4T": "pncA_GA-5A",
            "pncA_GTC-3GT": "pncA_GAC-5AC",
            "pncA_CGC442CGCGACGCGGTACGC": "pncA_CGC442CGCGACGCGGTACGC"
        },
    script:
        str(scripts_dir / "select_panel_variants.py")

# adjustments are explained in https://github.com/mbhall88/WHO-correspondence/issues/4#issuecomment-1096023422
rule construct_who_panel:
    input:
        vcf_dir=rules.extract_background_vcfs.output.vcf_dir,
        reference=h37rv,
        genbank=annotation,
        panel=rules.select_panel_variants.output.panel,
        resistance_json=rules.select_panel_variants.output.resistance_json,
    output:
        probes=panel_dir / "who2021/grading_{grade}/probes.fa",
    log:
        log_dir / "construct_who_panel/grading_{grade}.log",
    shadow:
        "shallow"
    container:
        containers["mykrobe"]
    params:
        db=".mongodb",
        db_name="who_panel",
        kmer_size=config["panel_kmer_size"],
    shell:
        """
        date +"[%Y-%m-%d %H:%M:%S] Starting the database..." > {log}
        mkdir -p {params.db} 2>> {log}
        mongod --quiet --dbpath {params.db} &> /dev/null &

        date +"[%Y-%m-%d %H:%M:%S] Adding background variants to database..." >> {log}
        for VCF in {input.vcf_dir}/*.vcf; do
          date +"[%Y-%m-%d %H:%M:%S] Adding variants for $VCF"
          mykrobe variants add -m cortex -f --db_name {params.db_name} "$VCF" {input.reference}
        done &>> {log}

        date +"[%Y-%m-%d %H:%M:%S] Making probes..." >> {log}
        mykrobe variants make-probes --db_name {params.db_name} -k {params.kmer_size} \
          -g {input.genbank} -t {input.panel} {input.reference} > {output.probes} 2>> {log}

        date +"[%Y-%m-%d %H:%M:%S] Shutting down the database..." >> {log}
        mongod --shutdown --dbpath {params.db} 2>> {log}

        date +"[%Y-%m-%d %H:%M:%S] Done!" >> {log}
        """
