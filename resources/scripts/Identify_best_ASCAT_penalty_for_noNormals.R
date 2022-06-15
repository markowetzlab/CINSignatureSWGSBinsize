# Find optimal ASCAT penalty based on ASCAT metadata when using ASCAT with matched normal as gold standard

library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(Cairo)

BASE="/Users/drews01/cnsigs2_revisions"
#META="/Users/drews01/cnsigs2_revisions/data/ASCAT_metadata_478_noNormals.tsv"
META="/Users/drews01/cnsigs2_revisions/data/ASCAT_metadata_478_WGS_downsampled_SNP6_QC.tsv"
#OUT="figures/478_vs_CELnoNorm_purity_ploidy_penalties.pdf"
OUT="figures/478_vs_WGStoSNP6withNorm_purity_ploidy_penalties.pdf"
GOLD=70
PENALTIES=c(35, 50, 70, 100, 140)
# tumour_mapd identical for all samples and penalties
# purity / ploidy and group investigated in first part
# segDiff not helpful as is a direct function of penalty
# LOH and frac_homo not helpful for deciding penalty as they give only information about what
# information we lose
# CRITERIA=c("tumour_mapd", "segDiff", "frac_homo", "purity", "ploidy", "LOH", "Group")
# TOL=0.1


meta = fread(META)

# Just by purity and ploidy - "group" column
lGroups = lapply(PENALTIES, function(pen) {
  
  # Could also be "group_penalty"
  dtOut = data.table(table(meta[[paste0("group_penalty", pen)]]))
  colnames(dtOut) = c("group", "N")
  dtOut$penalty = pen
  return(dtOut)
  
})

dtGroups = rbindlist(lGroups)
dtGroups$penalty = factor(dtGroups$penalty)
dtGroups$relative = dtGroups$N / 478

# ggplot(dtGroups, aes(x = group, y = relative, group = penalty, colour = penalty)) + geom_point()

# CEL files without normals
# p1 = ggplot(dtGroups, aes(x = penalty, y = relative, group = group, colour = group)) + 
#   geom_hline(yintercept = 0.192, colour = "grey80", linetype = "dashed") +
#   geom_hline(yintercept = 0.138, colour = "grey80", linetype = "dashed") +
#   geom_hline(yintercept = 0.086, colour = "grey80", linetype = "dashed") +
#   geom_hline(yintercept = 0.59, colour = "grey80", linetype = "dashed") +
#   geom_line() + geom_point() + labs(x = "ASCAT penalty", y = "Sample proportion") + 
#   scale_y_continuous(labels = scales::percent) + labs(colour = "Purity / ploidy classification")

# WGStoSNP6 with normals
p1 = ggplot(dtGroups, aes(x = penalty, y = relative, group = group, colour = group)) + 
  geom_hline(yintercept = 0.449, colour = "grey80", linetype = "dashed") +
  geom_hline(yintercept = 0.311, colour = "grey80", linetype = "dashed") +
  geom_hline(yintercept = 0.082, colour = "grey80", linetype = "dashed") +
  geom_hline(yintercept = 0.178, colour = "grey80", linetype = "dashed") +
  geom_line() + geom_point() + labs(x = "ASCAT penalty", y = "Sample proportion") + 
  scale_y_continuous(labels = scales::percent) + labs(colour = "Purity / ploidy classification")


cairo_pdf(file.path(BASE, OUT), height = 3.5, width = 5.5)
p1; dev.off()


# # Compare other criteria by a margin of 10% from gold standard
# PENALTIES = PENALTIES[ PENALTIES != GOLD ]
# lapply(CRITERIA, function(crit) {
#   
#   vGold = meta[[paste0(crit, "_penalty", GOLD)]]
#   lPens = lapply(PENALTIES, function(pen) {
#     
#     vPen = meta[[paste0(crit, "_penalty", pen)]]
#     
#   })
#   
# })


