library(dplyr)
library(methods)
library(GSVA)
library(GSEABase)
library(largeList)
library(genefilter)
library(stringi)
library(limma)
###For GSE13355
genesets <- getGmt("./GSVA/coreNetwork.gmt")

meta<-pdata %>% filter(characteristics_ch1 %in% c("involved skin from cases","uninvolved skin from cases"))
gene = as.matrix(exprSet4GSEA[,rownames(meta)])
meta$characteristics_ch1

nGrp1 <- 58   ###non_lesional
nGrp2 <- 58   ###lesional 
design <- cbind(sampleGroup1=c(rep(1, nGrp1), rep(0, nGrp2)), sampleGroup2=c(rep(0, nGrp1), rep(1, nGrp2)))
colnames(design) <- c("NL","L")

topMatrixGSVA <- gsva(gene, genesets, min.sz=10, max.sz=Inf, tau=1, method="gsva", abs.ranking=FALSE, verbose=TRUE, parallel.sz=12) 
adjPvalueCutoff <- 0.05   
fit <- lmFit(topMatrixGSVA, design)  
fit <- eBayes(fit)      
res <- decideTests(fit, p.value=adjPvalueCutoff)
summary(res)
write.csv(topMatrixGSVA, file = "GSE13355_GSVAScore.csv")

Score<-as.data.frame(t(topMatrixGSVA))
class(Score$coreNetwork) <- "numeric"
Score$group=c(rep("NL", nGrp1), rep("L", nGrp2))
Score$group<-factor(Score$group,levels = c("NL","L"))
library(ggpubr)
p1=ggplot(Score, aes(x=group, y=coreNetwork, color=group)) +
  geom_violin(width = 1) +
  geom_boxplot(width=0.1,position=position_dodge(0.9)) + scale_color_manual(values=c("#00468BFF", "#ED0000FF")) + theme_classic()+
  stat_compare_means(aes(label = ..p.signif..), method='t.test',paired=F, label.y.npc ="top",label.x.npc = 0.3,size=5, comparisons = list( c("NL", "L")))+
  labs(x = " ", y = "Core nework score", title = "Gudjonsson et al. 2009")
p1 + theme(axis.ticks = element_line(size = 1),
           axis.title = element_text(size = 23), 
           axis.text = element_text(size = 23), 
           axis.text.x = element_text(size = 23), 
           axis.text.y = element_text(size = 20)) +labs(x = NULL, colour = NULL)
ggsave("GSE13355_GSVAScore.pdf",width = 5, height = 5)

#Limma
