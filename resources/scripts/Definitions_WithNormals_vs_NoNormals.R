# Compare activities for the same samples from different high-throughput technologies

rm(list=ls(all=TRUE))

library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lsa)
library(patchwork)
theme_set(theme_tufte(base_size = 12, base_family = "Arial"))

# OUT1="~/cnsigs2_revisions/figures/4_Definitions_cel478_vs_cel478NoNorm_K10_ASCATcl.pdf"
OUT1="~/cnsigs2_revisions/figures/4_Definitions_cel478_vs_wgssnp6_K10.pdf"


mWithNorm = readRDS("~/cnsigs2_revisions/CINsignatures/Pancancer_signatures_Oct2020.rds")
# dtWithNorm = data.table(melt(mWithNorm))

# Add input from all others penalties

## ASCAT.cl WGS downsamples to SNP6 positions
vFiles = c("penalty50" = "~/cnsigs2_revisions/out/6_Signatures_E8610I.txt",
           "penalty70" = "~/cnsigs2_revisions/out/6_Signatures_G25CO5.txt",
           "penalty100" = "~/cnsigs2_revisions/out/6_Signatures_GTGB8P.txt",
           "penalty140" = "~/cnsigs2_revisions/out/6_Signatures_FAQTWZ.txt")
# 35 had no solution with K=10:
# "penalty35" = "~/cnsigs2_revisions/out/6_Signatures_VL8JB8.txt",

# ## ASCAT.sccl CELnoNorm
# vFiles = c("penalty35" = "~/cnsigs2_revisions/out/6_Signatures_MT5S79.txt",
#            "penalty50" = "~/cnsigs2_revisions/out/6_Signatures_JL88LL.txt",
#            "penalty70" = "~/cnsigs2_revisions/out/6_Signatures_TGTSA4.txt",
#            "penalty100" = "~/cnsigs2_revisions/out/6_Signatures_3JWBI9.txt",
#            "penalty140" = "~/cnsigs2_revisions/out/6_Signatures_1ZDVSR.txt")

# # ASCAT.cl CELnoNorm
# vFiles = c("penalty35" = "~/cnsigs2_revisions/out/6_Signatures_K0KUA5.txt",
#            "penalty50" = "~/cnsigs2_revisions/out/6_Signatures_A2H2K6.txt",
#            "penalty70" = "~/cnsigs2_revisions/out/6_Signatures_5N498Y.txt",
#            "penalty100" = "~/cnsigs2_revisions/out/6_Signatures_K6DGUK.txt",
#            "penalty140" = "~/cnsigs2_revisions/out/6_Signatures_HJUVT9.txt")

lNoNorm = list()
for(i in 1:length(vFiles)) {

    theseSigs = fread(vFiles[i])
    vFeats = theseSigs$V1
    mSigs = t(as.matrix(theseSigs[,-1]))
    colnames(mSigs) = vFeats
    rownames(mSigs) = paste0(rownames(mSigs), "_", names(vFiles[i]))
    lNoNorm[[length(lNoNorm)+1]] = mSigs
}
names(lNoNorm) = names(vFiles)

 
# mNoNorm = lNoNorm[["penalty70"]]
# dtNoNorm = data.table(melt(mNoNorm))

# ## Compare absolute values of signatures.
# allSigs = levels(dtWithNorm$Var2)
# lComp = lapply(allSigs, function(thisSig) {
#   
#   thisOld = dtWithNorm[ dtWithNorm$Var2 == thisSig, ]
#   thisNew = dtNoNorm[ dtNoNorm$Var2 == thisSig, ]
#   thisOld$valueNew = thisNew$value[ match(thisOld$Var1, thisNew$Var1) ]
#   return(thisOld)
#   
# })
# 
# dtComp = rbindlist(lComp)
# colnames(dtComp) = c("Sample", "Sig", "WithNormal", "NoNormal")
# dtComp$diff = dtComp$WithNormal - dtComp$NoNormal
# 
# # Todo: Mark purity, ploidy outliers
# p1 = ggplot(dtComp, aes(x = Sig, y = diff, group = Sample)) + 
#   geom_hline(yintercept = 0, colour = 'black') +
#   geom_hline(yintercept = 0.05, colour = 'red', linetype = "dashed") +
#   geom_hline(yintercept = -0.05, colour = 'red', linetype = "dashed") +
#   geom_jitter(height = 0, width = 0.3, alpha = 0.6) + 
#   ylab("WithNormal - NoNormal") + xlab("Signature") + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
#   
# cairo_pdf(OUT1, height = 4, width = 5); p1; dev.off()


### Compare cosine similarities
lCos = list()
for(i in 1:length(lNoNorm)) {
  
  thisNoNorm = lNoNorm[[i]]
  mWithNorm = mWithNorm[ ,match(colnames(thisNoNorm), colnames(mWithNorm)) ]
  
  vCos = sapply(1:nrow(mWithNorm), function(i) { 
    max(cosine(mWithNorm[i,], t(thisNoNorm)))
  })
  dtCos = data.table("sig" = rownames(mWithNorm), "maxcos" = vCos, "penalty" = names(lNoNorm[i]))
  lCos[[length(lCos) + 1]] = dtCos
  
}

dtCos = rbindlist(lCos)
dtCos$penalty = factor(dtCos$penalty, 
                  levels = c("penalty35", "penalty50", "penalty70", "penalty100", "penalty140"), 
                  labels = c(35, 50, 70, 100, 140))
dtCos$sig = factor(dtCos$sig, levels = rev(unique(dtCos$sig)))

# Categorise max cosine similarity
dtCos$category = "<0.5"
dtCos$category[ dtCos$maxcos >= 0.5 & dtCos$maxcos < 0.8 ] = "<0.8"
dtCos$category[ dtCos$maxcos >= 0.8 & dtCos$maxcos < 0.95 ] = "<0.95"
# dtCos$category[ dtCos$maxcos >= 0.9 & dtCos$maxcos < 0.95 ] = "<0.95"
dtCos$category[ dtCos$maxcos >= 0.95 ] = ">0.95"
dtCos$category = factor(dtCos$category, levels = c("<0.5", "<0.8", "<0.95", ">0.95"))

p1 = ggplot(dtCos, aes(x = penalty, y = sig, fill = category)) + geom_tile() +
  scale_x_discrete(drop=FALSE) + xlab("ASCAT Penalty") + ylab("Gold standard pan-cancer sigs") + 
  labs(fill = "Max cosine\nsimilarity") +
  # scale_fill_gradient2(low = "white", mid = "white", midpoint = 0.5, high = "#cb181d", limits = c(0,1))
  scale_fill_manual(values = c("#ffffcc", "#fed976", "#fd8d3c", "#e31a1c"), drop = FALSE)

cairo_pdf(OUT1, height = 4, width = 5); p1; dev.off()

