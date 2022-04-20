

## Library loading and data preparation


library(WGCNA)
load("wgcnaLesion.RData")
options(stringsAsFactors = FALSE)
datTraits <- coldata
datExpr <- rnaExpr
datExpr[1:6, 1:6]


# Multi-thread
enableWGCNAThreads(nThreads = 20)

# Soft power estimation
powers = c(seq(4,10,by=1), seq(12,20, by=2))
sft = pickSoftThreshold(datExpr, 
                        powerVector = powers, 
                        verbose = 2,
                        networkType = "signed")
powerTable <- sft[[2]]

pdf("softpower.pdf",width = 7,height = 5)
par(mfrow = c(1,2))
cex1 = 1
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
#abline(h=0.80,col="red")
# Mean connectivity ~ soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
# Get the scale-free network 
net = blockwiseModules(
    datExpr, 
    maxBlockSize = 20000, 
    corType = "bicor", # pearson|bicor
    power = 7, 
    networkType = "signed", # unsigned | signed | signed hybrid
    
    TOMType = "signed", # none | unsigned | signed
    saveTOMs = TRUE, 
    saveTOMFileBase = "blockwiseTOM", 
    minModuleSize = 30, 
    mergeCutHeight = 0.2, 
    numericLabels = TRUE, 
    nThreads = 0)

moduleColors = labels2colors(net$colors)
moduleLabels = net$colors
table(net$colors) 
table(moduleColors)


## Visualization


geneTree = net$dendrograms[[1]]
plotDendroAndColors(
    geneTree, 
    moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, 
    addGuide = TRUE)



## Clinical phenotyping data import and rownames check
sum(rownames(datTraits) != rownames(datExpr))

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# moduleEigengenes(MEs) calculation
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs)

# correlation between MEs and Phenotyping
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# heatmap plotting
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


plotMEpairs(MEs,y=datTraits$PASI, main = "PASI.Relationship between module eigengenes")
plotMEpairs(MEs,y=datTraits_swab$swab_total, main = "Local PASI.Relationship between module eigengenes")

## Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8,plotAdjacency = TRUE, signed = TRUE,
                      xLabelsAngle = 90)

## Plot the relationships among the eigengenes and PASI
PASI = as.data.frame(datTraits$PASI)
names(PASI) = "PASI"
MET = orderMEs(cbind(MEs, PASI))
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, plotAdjacency = TRUE, signed = TRUE,
                      xLabelsAngle = 90)

## Plot the relationships among the eigengenes and biopsy lesion scores
localPASI = as.data.frame(datTraits$swab_total) 
names(localPASI) = "local PASI"
MET = orderMEs(cbind(MEs, localPASI))
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8,plotAdjacency = TRUE, signed = TRUE,
                      xLabelsAngle = 90)


# Get the module and phenotyping names
modNames = substring(names(MEs), 3)
traitNames = names(datTraits)

# Get the MM(module membership matrix), correlation between genes and MEs
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p")) #Generalizing intramodular connectivity for all genes on the array

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                          nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# Get the GS (gene Significance), correlation between genes and phenotyping
geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                          nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="");
names(GSPvalue) = paste("p.GS.", traitNames, sep="");

geneInfo<-cbind(geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue)
write.table(geneInfo, file = "geneInfoLesion.txt", 
            sep = "\t", 
            quote = F)

save(powerTable, geneTree, MEs, moduleColors, moduleTraitCor,moduleTraitPvalue,moduleLabels, geneInfo, file = "NetworkConstruction-lesion.RData")


load("NetworkConstruction-lesion.RData")


