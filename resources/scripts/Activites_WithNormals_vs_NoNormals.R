# Compare activities for the same samples from different high-throughput technologies

library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lsa)
library(patchwork)
theme_set(theme_tufte(base_size = 12, base_family = "Arial"))

## Output files
# # No normal CEL files - newASCAT.sc
# OUT1="~/cnsigs2_revisions/figures/2_Activities_cel478_vs_cel478NoNorm_penalty70.pdf"
# OUT2="~/cnsigs2_revisions/figures/3_Activities_cel478_vs_cel478NoNorm_CosineSim.pdf"

# No normal CEL files - newASCAT
OUT1="~/cnsigs2_revisions/figures/2_Activities_cel478_vs_cel478NoNorm_penalty70_newASCATvsnewASCATsc.pdf"
OUT2="~/cnsigs2_revisions/figures/3_Activities_cel478_vs_cel478NoNorm_CosineSim_newASCATvsnewASCATsc.pdf"
OUT3="~/cnsigs2_revisions/figures/2_Activities_cel478_vs_cel478NoNorm_penalty70_ASCATcl.pdf"

# # WGS down SNP6 with normals
# OUT1="~/cnsigs2_revisions/figures/2_Activities_cel478_vs_wgs478withNorm_penalty70.pdf"
# OUT2="~/cnsigs2_revisions/figures/3_Activities_cel478_vs_wgs478withNorm_CosineSim.pdf"

# # ASCATsc Methylation
# OUT1="~/cnsigs2_revisions/figures/2_Activities_cel478_vs_Methylation.pdf"
# OUT2="~/cnsigs2_revisions/figures/3_Activities_cel478_vs_Methylation_CosineSim.pdf"


## Load gold standard (comparison signature activities)
mWithNorm = readRDS("~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures.rds")
dtWithNorm = data.table(melt(mWithNorm))

#### Add input from all others penalties
# ## CEL noNormals - oldASCAT vs newASCAT.sc
# vFiles = c("~/cnsigs2_revisions/out/NormalASCAT/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_noNormals_penalty35.rds",
#            "~/cnsigs2_revisions/out/NormalASCAT/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_noNormals_penalty50.rds",
#            "~/cnsigs2_revisions/out/NormalASCAT/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_noNormals_penalty70.rds",
#            "~/cnsigs2_revisions/out/NormalASCAT/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_noNormals_penalty100.rds",
#            "~/cnsigs2_revisions/out/NormalASCAT/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_noNormals_penalty140.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_CEL_noNormals_newASCAT_penalty35.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_CEL_noNormals_newASCAT_penalty50.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_CEL_noNormals_newASCAT_penalty70.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_CEL_noNormals_newASCAT_penalty100.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_CEL_noNormals_newASCAT_penalty140.rds"
#            )

# ## CEL noNormals - oldASCAT vs newASCAT
# vFiles = c("~/cnsigs2_revisions/out/NormalASCAT/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_noNormals_penalty35.rds",
#            "~/cnsigs2_revisions/out/NormalASCAT/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_noNormals_penalty50.rds",
#            "~/cnsigs2_revisions/out/NormalASCAT/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_noNormals_penalty70.rds",
#            "~/cnsigs2_revisions/out/NormalASCAT/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_noNormals_penalty100.rds",
#            "~/cnsigs2_revisions/out/NormalASCAT/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_noNormals_penalty140.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty35.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty50.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty70.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty100.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty140.rds"
#            )

## CEL noNormals - ASCAT.cl
vFiles = c("~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty35.rds",
           "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty50.rds",
           "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty70.rds",
           "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty100.rds",
           "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty140.rds")
           

# ## CEL noNormals - newASCATsc vs newASCAT
# vFiles = c("~/cnsigs2_revisions/out/NewASCATsc_CELnoNorm/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_CEL_noNormals_newASCAT_penalty35.rds",
#            "~/cnsigs2_revisions/out/NewASCATsc_CELnoNorm/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_CEL_noNormals_newASCAT_penalty50.rds",
#            "~/cnsigs2_revisions/out/NewASCATsc_CELnoNorm/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_CEL_noNormals_newASCAT_penalty70.rds",
#            "~/cnsigs2_revisions/out/NewASCATsc_CELnoNorm/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_CEL_noNormals_newASCAT_penalty100.rds",
#            "~/cnsigs2_revisions/out/NewASCATsc_CELnoNorm/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_CEL_noNormals_newASCAT_penalty140.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty35.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty50.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty70.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty100.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty140.rds"
# )


# ## WGSdownSNP6 withNormals
# vFiles = c("~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGSdownSNP6_withNormals_penalty35.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGSdownSNP6_withNormals_penalty50.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGSdownSNP6_withNormals_penalty70.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGSdownSNP6_withNormals_penalty100.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGSdownSNP6_withNormals_penalty140.rds")

# ## PCAWG Methylation
# vFiles = c("~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_CEL_noNormals_newASCAT_penalty50.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_CEL_noNormals_newASCAT_penalty70.rds",
#            "~/cnsigs2_revisions/out/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_PCAWG_TCGA_ASCATsc_Methylation.rds")

## Load files
lNoNorm = lapply(vFiles, function(thisFile) readRDS(thisFile))

# # CEL files noNorm with comparison to other ASCAT versions
# names(lNoNorm) = c("penalty35NewASCATsc", "penalty50NewASCATsc", "penalty70NewASCATsc", 
#                    "penalty100NewASCATsc", "penalty140NewASCATsc",
#                    "penalty35NewASCAT", "penalty50NewASCAT", "penalty70NewASCAT",
#                    "penalty100NewASCAT", "penalty140NewASCAT")

# CEL files noNorm
names(lNoNorm) = c("penalty35", "penalty50", "penalty70", 
                   "penalty100", "penalty140")


# # Methylation
# names(lNoNorm) = c("penalty50", "penalty70", "Methylation")


## Plot 1: Compare activities per signature directly
# mNoNorm = lNoNorm[["penalty70NewASCATsc"]]
mNoNorm = lNoNorm[["penalty70"]]
dtNoNorm = data.table(melt(mNoNorm))

## Compare absolute values of signatures.
allSigs = levels(dtWithNorm$Var2)
lComp = lapply(allSigs, function(thisSig) {
  
  thisOld = dtWithNorm[ dtWithNorm$Var2 == thisSig, ]
  thisNew = dtNoNorm[ dtNoNorm$Var2 == thisSig, ]
  thisOld$valueNew = thisNew$value[ match(thisOld$Var1, thisNew$Var1) ]
  return(thisOld)
  
})

dtComp = rbindlist(lComp)
colnames(dtComp) = c("Sample", "Sig", "WithNormal", "NoNormal")
dtComp$diff = dtComp$WithNormal - dtComp$NoNormal

# Todo: Mark purity, ploidy outliers
p1 = ggplot(dtComp, aes(x = Sig, y = diff, group = Sample)) + 
  geom_hline(yintercept = 0, colour = 'black') +
  geom_hline(yintercept = 0.05, colour = 'red', linetype = "dashed") +
  geom_hline(yintercept = -0.05, colour = 'red', linetype = "dashed") +
  geom_jitter(height = 0, width = 0.3, alpha = 0.6) + 
  ylab("WithNormal - NoNormal") + xlab("Signature") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
cairo_pdf(OUT1, height = 4, width = 5); p1; dev.off()


### Compare cosine similarities
lCos = sapply(lNoNorm, function(thisNoNorm) {
  
  mWithMatched = mWithNorm[ match(rownames(thisNoNorm), rownames(mWithNorm)), ]
  vCos = sapply(1:nrow(mWithMatched), function(i) { cosine(mWithMatched[i,], thisNoNorm[i,]) })
  
  
}, simplify = FALSE, USE.NAMES = TRUE)

dtCos = data.table(melt(lCos))
# CEL files noNorm
dtCos$L1 = factor(dtCos$L1,
                  levels = c("penalty35", "penalty50", "penalty70",
                             "penalty100", "penalty140"),
                  labels = c("35", "50", "70", "100", "140"))

# # CEL files noNorm - incl. comparison
# dtCos$L1 = factor(dtCos$L1,
#                   levels = c("penalty35NewASCATsc", "penalty35NewASCAT", 
#                              "penalty50NewASCATsc", "penalty50NewASCAT", 
#                              "penalty70NewASCATsc", "penalty70NewASCAT", 
#                              "penalty100NewASCATsc", "penalty100NewASCAT",
#                              "penalty140NewASCATsc", "penalty140NewASCAT"),
#                   labels = c("35sc", "35ASCAT", "50sc", "50ASCAT", "70sc", "70ASCAT", "100sc", "100ASCAT", "140sc", "140ASCAT"))

# Methylation
# dtCos$L1 = factor(dtCos$L1, 
#                   levels = c("penalty50", "penalty70", "Methylation"), 
#                   labels = c(50, 70, "Methylation"))

p2a = ggplot(dtCos, aes(x = value, colour = L1)) + geom_density(alpha = 0.5, adjust = 0.5) + 
  xlab("Cosine similarity to matched normal") + ylab("Density") + labs(colour = "Penalty") + 
  ggtitle("Density adjust parameter = 0.5") + theme(legend.position = "none") #+
  # scale_colour_brewer(palette = "Paired")
p2b = ggplot(dtCos, aes(x = value, colour = L1)) + geom_density(alpha = 0.5, adjust = 1) + 
  xlab("Cosine similarity to matched normal") + ylab("Density") + labs(colour = "Penalty") + 
  ggtitle("Density adjust parameter = 1 (default)") + theme(legend.position = "none") #+
  # scale_colour_brewer(palette = "Paired")
p2c = ggplot(dtCos, aes(x = value, colour = L1)) + geom_density(alpha = 0.5, adjust = 2) + 
  xlab("Cosine similarity to matched normal") + ylab("Density") + labs(colour = "Penalty") + 
  ggtitle("Density adjust parameter = 2") #+
  # scale_colour_brewer(palette = "Paired")

p2 = p2a + p2b + p2c + plot_layout(widths = c(0.3, 0.3, 0.4))
cairo_pdf(OUT3, height = 3.5, width = 14); p2; dev.off()

