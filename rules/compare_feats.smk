rule compare_feats:
    input:
        feats=expand(OUT_DIR+"sWGS_binsize/{bin}kb/features_signatures/{bin}kb_features.rds",bin=BIN_VALS)
    output:
        OUT_DIR+"sWGS_binsize/seg_size_MannWhit_test.tsv"
    params:
        outdir=OUT_DIR
    script:
        "../scripts/compare_features.R"
