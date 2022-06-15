# Clean env
args = commandArgs(trailingOnly=TRUE)

#load libraries
library(QDNAseqmod)
library(Biobase)
library(ggplot2)
library(stringr)
suppressWarnings(library(doMC))
suppressWarnings(library(foreach))

#get segs table
getSegTable<-function(x)
{
  dat<-x
  sn<-Biobase::assayDataElement(dat,"segmented")
  fd <- Biobase::fData(dat)
  fd$use -> use
  fdfiltfull<-fd[use,]
  sn<-sn[use,]
  segTable<-c()
  for(c in unique(fdfiltfull$chromosome))
  {
    snfilt<-sn[fdfiltfull$chromosome==c]
    fdfilt<-fdfiltfull[fdfiltfull$chromosome==c,]
    sn.rle<-rle(snfilt)
    starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
    ends <- cumsum(sn.rle$lengths)
    lapply(1:length(sn.rle$lengths), function(s) {
      from <- fdfilt$start[starts[s]]
      to <- fdfilt$end[ends[s]]
      segValue <- sn.rle$value[s]
      c(fdfilt$chromosome[starts[s]], from, to, segValue)
    }) -> segtmp
    segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),stringsAsFactors=F)
    segTable<-rbind(segTable,segTableRaw)
  }
  segTable$sample <- rep(sampleNames(dat),times=nrow(segTable))
  colnames(segTable) <- c("chromosome", "start", "end", "segVal","sample")
  segTable
}

qc.data <- read.table(snakemake@input[["meta"]],header = T,sep = "\t")
output_dir <- snakemake@params[["outdir"]]
bin <- as.numeric(snakemake@params[["bin"]])
cores <- as.numeric(snakemake@threads)
registerDoMC(cores)

rds.filename <- snakemake@input[["rds"]]
rds.list <- lapply(rds.filename,FUN=function(x){readRDS(x)})

collapse_rds <- function(rds.list){
  comb <- rds.list[[1]][[1]]
  if(length(rds.list) > 1){
    for(i in 2:length(rds.list)){
      add <- rds.list[[i]][[1]]
      comb <- combine(comb,add)
    }
    rds.obj <- comb
  }
  return(rds.obj)
}

# Combine and load rds objects
rds.rel <- collapse_rds(rds.list)

if(!dir.exists(paste0(output_dir,"sWGS_binsize/",bin,"kb/abs_cn_rds/plots"))){
	dir.create(paste0(output_dir,"sWGS_binsize/",bin,"kb/abs_cn_rds/plots"))
}

# convert depth to abs cn
depthtocn<-function(x,purity,seqdepth) #converts readdepth to copy number given purity and single copy depth
{
  (x/seqdepth-2*(1-purity))/purity
}

# List samples
samples <- qc.data[which(qc.data$submitted_donor_id %in% colnames(rds.rel)),]

# Add pheno information
pData(rds.rel)$purity <- samples$purity[match(pData(rds.rel)$name,samples$submitted_donor_id)]
pData(rds.rel)$ploidy <- samples$ploidy[match(pData(rds.rel)$name,samples$submitted_donor_id)]

# Generate abs plot and table of fits
res <- data.frame(matrix(ncol = 9, nrow = 0))
abs_profiles <- rds.rel[fData(rds.rel)$use,]
# For each
for(sample in pData(rds.rel)$name){
  # Index and subselect sample
  ind <- which(colnames(rds.rel)==sample)
  relcn <- rds.rel[,ind]
  to_use <- fData(relcn)$use #
  relcn <- relcn[to_use,]
  smooth.bool <- FALSE
  # Extract cn and ploidy
  copynumber <- assayDataElement(relcn,"copynumber")
  rel_ploidy <- mean(copynumber,na.rm=T)
  ploidy <- pData(relcn)$ploidy
  purity <- pData(relcn)$purity
  cellploidy <- ploidy*purity+2*(1-purity)
  seqdepth <- rel_ploidy/cellploidy

  # Extract CN and Segs
  cn <- assayDataElement(relcn,"copynumber")
  seg <- assayDataElement(relcn,"segmented")
  
  # Convert to abs
  abs_cn <- depthtocn(cn,purity,seqdepth)
  abs_seg <- depthtocn(seg,purity,seqdepth)
  assayDataElement(relcn,"copynumber") <- abs_cn
  assayDataElement(relcn,"segmented") <- abs_seg
  # Add to abs RDS
  assayDataElement(abs_profiles,"copynumber")[,ind] <- abs_cn
  assayDataElement(abs_profiles,"segmented")[,ind] <- abs_seg
  # Add TP53 info
  TP53cn<-round(depthtocn(seg[73504],purity,seqdepth),1) # to 1 decimal place / altered to correct bin value
  expected_TP53_AF <- TP53cn*purity/(TP53cn*purity+2*(1-purity))
  TP53freq <-NA
  # Add patient-level info
  pat <- sample
  res <- rbind(res,matrix(c(sample,pat,ploidy,purity,TP53cn,round(expected_TP53_AF,2),TP53freq,NA,NA),nrow = 1,ncol = 9))
  
  # Y axis range
  if(ploidy>5){
    yrange=15
  } else {
    yrange=10
  }
  # Plot abs fit
  png(paste0(output_dir,"sWGS_binsize/",bin,"kb/abs_cn_rds/plots/",sample,".png"), w= 8, h = 6, unit="in", res = 250)
  par(mfrow = c(1,1))
  plot(relcn,doCalls=FALSE,showSD=TRUE,logTransform=FALSE,ylim=c(0,yrange),ylab="Absolute tumour CN",
       main=paste(sample, " eTP53=",round(expected_TP53_AF,2),
                  " AF=", round(TP53freq,2),
                  " p=",round(purity,2),
                  " pl=",round(ploidy,2),
                  sep=""))
  abline(h=1:9, col = "blue")
  dev.off()
}

segTables <- do.call(rbind,lapply(sampleNames(abs_profiles),FUN=function(y){
        getSegTable(abs_profiles[,y])
	})
)
write.table(segTables,paste0(output_dir,"sWGS_binsize/",bin,"kb/abs_cn_rds/",bin,"kb_abs_segTable.tsv"),sep = "\t",quote=F,row.names=FALSE)
saveRDS(segTables,paste0(output_dir,"sWGS_binsize/",bin,"kb/abs_cn_rds/",bin,"kb_abs_segTable.rds"))

# Annotated and rename table
colnames(res) <- c("SAMPLE_ID","PATIENT_ID","ploidy","purity","TP53cn","expected_TP53_AF","TP53freq","use","notes")
res <- data.frame(res,stringsAsFactors = F)

# Save rds
saveRDS(abs_profiles,file=paste0(output_dir,"sWGS_binsize/",bin,"kb/abs_cn_rds/",bin,"kb_ds_absCopyNumber.rds"))

#write table of fits
write.table(res,paste0(output_dir,"sWGS_binsize/",bin,"kb/abs_cn_rds/",bin,"kb_ds_abs_fits.tsv"),sep = "\t",quote=F,row.names=FALSE)

