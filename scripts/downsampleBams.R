args = commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(QDNAseqmod))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(dplyr))

bam_in <- snakemake@input[["bam"]]
meta <- snakemake@input[["meta"]]
outdir <- snakemake@params[["outdir"]]
bin <- as.numeric(snakemake@params[["bin"]])
outname <- snakemake@output[[1]]
sample_name <- snakemake@params[["sample"]]

sample_meta <- read.table(file = meta,header = T,sep = "\t",na.strings = "")

purity <- sample_meta$purity[sample_meta$submitted_donor_id == sample_name]
ploidy <- sample_meta$ploidy[sample_meta$submitted_donor_id == sample_name]

bins<-getBinAnnotations(binSize = bin)
nbins_ref_genome <- nrow(bins@data[bins@data$use == TRUE,])

downsample_RD <- (((2*(1-purity)+purity*ploidy)/(ploidy*purity))/purity)*15*(2*(1-purity)+purity*ploidy)*nbins_ref_genome*(1/0.91)

cmd.fullRD <- paste("samtools view -c -F 260 ",bam_in)
full_rd <- as.numeric(system(cmd.fullRD,intern=T))

perc <- downsample_RD / full_rd
print(perc)

cmd.downsample <- paste("samtools view -s ", perc," -b ",bam_in," > ",outname)
cmd.index <- paste0("samtools index ",outname)
   
system(cmd.downsample)
system(cmd.index)
