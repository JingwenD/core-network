library(RegEnrich)
load("~/Project/2_RNAseqSkin/RegEnrich/BLN4RegEnrich.RData")
data(TFs)
head(TFs)
colData$Lesion <-as.factor(colData$Lesion)
colData$Lesion <-relevel(colData$Lesion,"NL")

object = RegenrichSet(expr = BLN,
                      colData = colData,
                      method = "LRT_DESeq2", 
                      minMeanExpr = 0.01,
                      design = ~ Lesion,
                      reduced = ~1,
                      networkConstruction = "GRN")
print(object)

## ----regenrich_diffExpr-------------------------------------------------------------------------------------
object = regenrich_diffExpr(object)
print(object)
print(results_DEA(object))

library(BiocParallel)
# on non-Windows computers (use 2 workers)
bpparam = register(MulticoreParam(workers = 6, RNGseed = 1234))
# on Windows computers (use 2 workers)
#bpparam = register(SnowParam(workers = 2, RNGseed = 1234))

object_grn = regenrich_network(object, networkConstruction = "GRN",
                               BPPARAM = bpparam, minR = 0.2)
print(object_grn)
print(results_topNet(object_grn))
save(object_grn, file = "grn_BLN.Rdata")

## ----regenrich_enrich_GSEA----------------------------------------------------------------------------------
load(file = "grn_BLN.Rdata")
set.seed(123)
object_grn_GSEA = regenrich_enrich(object_grn, enrichTest = "GSEA", nperm = 10000)
print(results_enrich(object_grn_GSEA))
grn_enrich_GSEA = results_enrich(object_grn_GSEA)@allResult
grn_enrich_GSEA = slot(results_enrich(object_grn_GSEA), "allResult")
head(grn_enrich_GSEA[, 1:6])



## ----regenrich_rankScore------------------------------------------------------------------------------------
grn_FET_rankScore = regenrich_rankScore(object_grn_FET)
results_score(grn_FET_rankScore)

## ----regenrich_rankScore------------------------------------------------------------------------------------
grn_GSEA_rankScore = regenrich_rankScore(object_grn_GSEA)
results_score(grn_GSEA_rankScore)

write.csv(grn_GSEA_rankScore@resScore,file = "grn_GSEA_rankScore.csv")
