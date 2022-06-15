# Find optimal ASCAT penalty based on ASCAT metadata when using ASCAT with matched normal as gold standard

# This script is just an extention for CEL files without normals based on Maxime's work.

## Functions
addGroupColumn = function(dtMeta, purDiff = 0.05, ploiDiff = 0.1){
  
  # First purity
  dtMeta$groupPur = TRUE
  dtMeta$groupPur[ abs( dtMeta$oirinalPurity - dtMeta$purity ) > purDiff ] = FALSE
  
  # Second ploidy
  dtMeta$groupPloid = TRUE
  dtMeta$groupPloid[ abs( dtMeta$originalTumourPloidy - dtMeta$tumour_ploidy ) > ploiDiff ] = FALSE
  
  # Make group column
  dtMeta$group = "Different purity and different ploidy"
  dtMeta$group[ dtMeta$groupPur & dtMeta$groupPloid ] = "Similar purity and similar ploidy "
  dtMeta$group[ ! dtMeta$groupPur & dtMeta$groupPloid ] = "Different purity and similar ploidy "
  dtMeta$group[ dtMeta$groupPur & ! dtMeta$groupPloid ] = "Similar purity and different ploidy "
  
  # Clean up (also Maxime's columns for easier overview)
  dtMeta[ , c("groupPloid", "groupPur", "isSamePur", "isSamePsi", 
              "isSamePurRelaxed", "isSamePsiRelaxed" ) ] = NULL
  
  return(dtMeta)
  
}

## Boilerplate
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(Cairo)

BASE="/Users/drews01/cnsigs2_revisions"
OUT=file.path(BASE,"figures/478_vs_CELnoNorm_purity_ploidy_penalties_newASCATsc.pdf")
vFiles=c(file.path(BASE, "data/NewASCATsc_metadata_CELnoNorm/summary_purity_ploidy_penalty35.txt"),
         file.path(BASE, "data/NewASCATsc_metadata_CELnoNorm/summary_purity_ploidy_penalty50.txt"),
         file.path(BASE, "data/NewASCATsc_metadata_CELnoNorm/summary_purity_ploidy_penalty70.txt"),
         file.path(BASE, "data/NewASCATsc_metadata_CELnoNorm/summary_purity_ploidy_penalty100.txt"),
         file.path(BASE, "data/NewASCATsc_metadata_CELnoNorm/summary_purity_ploidy_penalty140.txt"))


# Load data and identify purity / ploidy groups. Also identify all sample names.
lGroups = lapply(vFiles, function(thisFile) {
  
  dtMeta = fread(thisFile)
  pen = strsplit(strsplit(basename(thisFile), '_')[[1]][4], "\\.")[[1]][1]
  dtMeta = addGroupColumn(dtMeta)
  
  dtOut = data.table(melt(table(dtMeta$group)))
  colnames(dtOut) = c("group", "N")
  dtOut$penalty = pen
  return(dtOut)
  
})

dtGroups = rbindlist(lGroups)
dtGroups$penalty = factor(substr(dtGroups$penalty, 8, 11), levels = c(35, 50, 70, 100, 140))
dtGroups$relative = dtGroups$N / 478

# Plot
p1 = ggplot(dtGroups, aes(x = penalty, y = relative, group = group, colour = group)) + 
  geom_hline(yintercept = 0.567, colour = "grey80", linetype = "dashed") +
  geom_hline(yintercept = 0.22, colour = "grey80", linetype = "dashed") +
  geom_hline(yintercept = 0.138, colour = "grey80", linetype = "dashed") +
  geom_hline(yintercept = 0.077, colour = "grey80", linetype = "dashed") +
  geom_line() + geom_point() + labs(x = "ASCAT penalty", y = "Sample proportion") + 
  scale_y_continuous(labels = scales::percent) + labs(colour = "Purity / ploidy classification")

cairo_pdf(OUT, height = 3.5, width = 5.5)
p1; dev.off()


