# Merge TCGA segment files to one summary file

library(data.table)

FOLDER="~/cnsigs2_revisions/data/NewASCAT_CELnoNorm/profiles_ASCAT_70"
OUTPUT="~/cnsigs2_revisions/data/NewASCAT_CELnoNorm/NewASCAT_478TCGA_CELnoNorm_penalty70.rds"
MAXIME=FALSE
allFiles = list.files(FOLDER, full.names = TRUE)

lOut = lapply(allFiles, function(thisFile) {

  if(file.size(thisFile)==1) { return(NULL) }
  raw = fread(thisFile)
  if(MAXIME) {
    # Remove segments with NA
    raw = raw[ ! is.na(raw$total_copy_number), ]
    
    # Litte bug fix for now
    raw$chromosome[ is.na(raw$chromosome) ] = "X"
    
    sample = substr(basename(thisFile), 1, 12)
    out = raw[, c("chromosome", "start", "end", "total_copy_number_logr")]
    out$sample = sample
    colnames(out) = c("chromosome", "start", "end", "segVal", "sample")
    return(out)
  } else {
    raw$segVal = raw$nAraw + raw$nBraw
    out = raw[, c("chr", "startpos", "endpos", "segVal", "sample")]
    colnames(out) = c("chromosome", "start", "end", "segVal", "sample")
    return(out)
    
  }
})

dtOut = rbindlist(lOut)

# Replace "." in sample names with "-"
if( sum(grepl(".", dtOut$sample)) > 0) {
  dtOut$sample = gsub(pattern = "\\.", replacement = "-", dtOut$sample)
}
  
saveRDS(dtOut, OUTPUT)
