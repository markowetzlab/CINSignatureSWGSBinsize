include: "rules/common.smk"

def get_bam(wildcards):
    files = list(samplesheet.loc[(wildcards.sample), ["file"]])
    return files

SAMPLES = list(samplesheet['submitted_donor_id'].unique())

localrules: all

rule all:
    input:
        OUT_DIR+"sWGS_binsize/seg_size_MannWhit_test.tsv"

include: "rules/symlink.smk"
include: "rules/downsample.smk"
include: "rules/rel_rds.smk"
include: "rules/rel_to_abs.smk"
include: "rules/gen_feat_sigs.smk"
include: "rules/compare_feats.smk"
