rule sym_link:
    input:
        bam=get_bam
    output:
        expand(OUT_DIR+"sWGS_binsize/{{bin}}kb/bams/{{sample}}.bam")
    threads: 1
    shell:
        "ln -s {input.bam} {output}"
