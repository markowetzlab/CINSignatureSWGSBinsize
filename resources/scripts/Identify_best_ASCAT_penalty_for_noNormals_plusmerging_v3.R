# Find optimal ASCAT penalty based on ASCAT metadata when using ASCAT with matched normal as gold standard

# This script is just an extention for CEL files without normals based on Maxime's work.

## Boilerplate
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(Cairo)

BASE="/Users/drews01/cnsigs2_revisions"
OUT=file.path(BASE,"figures/478_vs_CELnoNorm_purity_ploidy_penalties_newASCAT.pdf")
META=file.path(BASE, "data/ASCAT_metadata_TCGA_PCAWG_with_dCIN.txt")
P70=file.path(BASE, "data/NewASCAT_metadata_CELnoNorm/summary_purity_ploidy_penalty70_ASCAT.txt")
POTHER=file.path(BASE, "data/NewASCAT_metadata_CELnoNorm/summary_purity_ploidy_penalties_ASCAT.txt")


# Merge metadata of new ASCAT runs
dtP70 = fread(P70)
dtP70$penalty = 70
dtOther = fread(POTHER)
dtOther$samplename = substr(dtOther$samplename, 1, 12)
dtNew = rbind(dtP70, dtOther)

# Load original metadata
meta = fread(META)
dtMeta = meta[, c("patient", "purity", "ploidy") ]

# Transfer original purity / ploidy
dtNew$oriPurity = dtMeta$purity[ match(dtNew$samplename, dtMeta$patient) ]
dtNew$oriPloidy = dtMeta$ploidy[ match(dtNew$samplename, dtMeta$patient) ]

# Assign groups based on Tom's stricter criteria
purDiff = 0.05
ploiDiff = 0.1
# First purity
dtNew$groupPur = TRUE
dtNew$groupPur[ abs( dtNew$oriPurity - dtNew$purity ) > purDiff ] = FALSE

# Second ploidy
dtNew$groupPloid = TRUE
dtNew$groupPloid[ abs( dtNew$oriPloidy - dtNew$ploidy ) > ploiDiff ] = FALSE

# Make group column
dtNew$group = "Different purity and different ploidy"
dtNew$group[ dtNew$groupPur & dtNew$groupPloid ] = "Similar purity and similar ploidy "
dtNew$group[ ! dtNew$groupPur & dtNew$groupPloid ] = "Different purity and similar ploidy "
dtNew$group[ dtNew$groupPur & ! dtNew$groupPloid ] = "Similar purity and different ploidy "

# Summarise groups per penalty
dtOut = data.table(melt(table(dtNew$group, dtNew$penalty)))
dtSum = aggregate(value ~ Var2, dtOut, sum)
dtOut$all = dtSum$value[ match(dtOut$Var2, dtSum$Var2) ]
dtOut$prop = dtOut$value / dtOut$all

# Plot results
dtOut$Var2 = factor(dtOut$Var2)
p1 = ggplot(dtOut, aes(x = Var2, y = prop, group = Var1, colour = Var1)) + 
  geom_hline(yintercept = 0.527, colour = "grey80", linetype = "dashed") +
  geom_hline(yintercept = 0.227, colour = "grey80", linetype = "dashed") +
  geom_hline(yintercept = 0.087, colour = "grey80", linetype = "dashed") +
  geom_hline(yintercept = 0.173, colour = "grey80", linetype = "dashed") +
  geom_line() + geom_point() + labs(x = "ASCAT penalty", y = "Sample proportion") + 
  scale_y_continuous(labels = scales::percent) + labs(colour = "Purity / ploidy classification")

cairo_pdf(OUT, height = 3.5, width = 5.5)
p1; dev.off()
