rule relRDS:
    input:
        bam=expand(OUT_DIR+"sWGS_binsize/{{bin}}kb/downsampled_bams/{{sample}}.bam"),
        meta=config["samples"]
    output:
        OUT_DIR+"sWGS_binsize/{bin}kb/relative_cn_rds/{sample}_{bin}kb_relSmoothedCN.rds"
    params:
        outdir=OUT_DIR,
        bin="{bin}",
        sample="{sample}"
    script:
        "../scripts/qdnaseq_mod_ds.R"
