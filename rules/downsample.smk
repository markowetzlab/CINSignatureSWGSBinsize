rule downsample:
    input:
        bam=OUT_DIR+"sWGS_binsize/{bin}kb/bams/{sample}.bam",
        meta=config["samples"]
    output:
        OUT_DIR+"sWGS_binsize/{bin}kb/downsampled_bams/{sample}.bam"
    params:
        bin="{bin}",
        sample="{sample}",
        outdir=OUT_DIR
    script:
        "../scripts/downsampleBams.R"
