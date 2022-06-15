rule feat_sigs:
    input:
        OUT_DIR+"sWGS_binsize/{bin}kb/abs_cn_rds/{bin}kb_abs_segTable.rds"
    output:
        OUT_DIR+"sWGS_binsize/{bin}kb/features_signatures/{bin}kb_features.rds"
    params:
        bin="{bin}",
        outdir=OUT_DIR
    script:
        "../scripts/gen_sigs_feats.R"
