rule rel_to_abs:
    input:
        rds=expand(OUT_DIR+"sWGS_binsize/{{bin}}kb/relative_cn_rds/{sample}_{{bin}}kb_relSmoothedCN.rds",sample=SAMPLES),
        meta=config["samples"]
    output:
        tsv=OUT_DIR+"sWGS_binsize/{bin}kb/abs_cn_rds/{bin}kb_ds_abs_fits.tsv",
        rds=OUT_DIR+"sWGS_binsize/{bin}kb/abs_cn_rds/{bin}kb_abs_segTable.rds"
    params:
        outdir=OUT_DIR,
        bin="{bin}"
    script:
        "../scripts/qdnaseq_rel_to_abs.R"

