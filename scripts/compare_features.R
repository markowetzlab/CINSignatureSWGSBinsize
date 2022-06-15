## Compare cn feats

library(ggplot2)
library(ggridges)
library(stringr)

## define extracted feature distributions
input_feats <- snakemake@input[["feats"]]
ref_feats <- readRDS("resources/CINsignatures/2_TCGA_PCAWG_ECNF.rds")
OUT_DIR <- snakemake@params[["outdir"]]
OUT_DIR <- paste0(OUT_DIR,"sWGS_binsize/")

input_names <- gsub(basename(input_feats),pattern = "_features.rds",replacement = "")

## Load extracted feature distributions
feat_table <- do.call(rbind,lapply(input_feats, FUN = function(x){
  input_name <- gsub(basename(x),pattern = "_features.rds",replacement = "")
  rds <- readRDS(x)
  rds.names <- names(rds)
  rds_sub <- do.call(rbind,lapply(rds.names, FUN = function(y){
      tab <- as.data.frame(rds[[y]])
      tab$feat <- rep(y,times=nrow(tab))
      colnames(tab) <- c("ID","value","feat")
      return(tab)
    }))
  rds_sub$bin <- rep(input_name,times=nrow(rds_sub))
  return(rds_sub)
  }))

## Subselect samples from total set
input_samples <- unique(feat_table$ID)

ref.names <- names(ref_feats)
ref_feats <- do.call(rbind,lapply(ref.names, FUN = function(y){
  tab <- as.data.frame(ref_feats[[y]])
  tab <- tab[tab$ID %in% input_samples,]
  tab$feat <- rep(y,times=nrow(tab))
  colnames(tab) <- c("ID","value","feat")
  return(tab)
}))
ref_feats$bin <- rep("ref",times=nrow(ref_feats))

feat_data <- rbind(feat_table,ref_feats)
feat_data$bin <- factor(feat_data$bin,levels = str_sort(numeric = T,unique(feat_data$bin)))

feat_data <- feat_data[feat_data$value < 10000000,]

png(paste0(OUT_DIR,"seg_size_plot_full.png"),width = 10,height = 8,units = "in",res = 600)
ggplot(feat_data[feat_data$feat == "segsize",]) +
  geom_density_ridges2(aes(x = value,y = bin,fill=bin),alpha=0.5) +
  geom_vline(data=feat_data[feat_data$feat == "segsize" & feat_data$bin == "ref",],aes(xintercept=median(value))) +
  labs(title = "Segment size distribution") +
  theme_ridges()
dev.off()

get_seg_size_signif <- function(data=NULL){
  tab <- c()
  s <- "segsize"
  ref.data <- data$value[data$bin == "ref" & data$feat == s]
  for(i in input_names){
        test.data <- data$value[data$bin == i & data$feat == s]
        p <- wilcox.test(x = ref.data,y = test.data)
        tabrow <- c(feat=s,ref="ref",test=i,p.value=p$p.value)
        tab <- rbind(tab,tabrow)
  }
  return(as.data.frame(tab,row.names = NULL))
}

print(feat_data)
mann_whit_tab <- get_seg_size_signif(data = feat_data)
write.table(mann_whit_tab,file = paste0(OUT_DIR,"seg_size_MannWhit_test.tsv"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
