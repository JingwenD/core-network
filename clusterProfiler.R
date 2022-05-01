library(clusterProfiler)
library(org.Hs.eg.db)
load("DEA.RData")
GO05 <- enrichGO(df_05$gene_name, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.2, keyType = 'SYMBOL') 
GO05@result<-GO05@result[order(GO05@result$Count,decreasing = T),]
pdf("GO_LvsNL_1.pdf", width=7, height=3.5)
dotplot(GO05,showCategory=10)
dev.off()
pdf(file="GOAll_LvsNL.svg",width=9,height=8)
dotplot(GO05, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
dev.off()
write.csv(GO05, file = "GO_LvsNL.csv")
EnZid05 = bitr(rownames(df_05),                #Gene name convert
               fromType="ENSEMBL", #input SYMBOL format, can be "ENSEMBL" or "SYMBOL"
               toType="ENTREZID",  # output ENTERZID format
               OrgDb="org.Hs.eg.db")

require(ReactomePA)
library("stringr")
RAT05 <- enrichPathway(gene=EnZid05$ENTREZID, pvalueCutoff=0.05, readable=T,organism = "human")
write.csv(RAT05,file="Reactome05.csv",quote=F)
RAT05@result<-RAT05@result[order(RAT05@result$Count,decreasing = T),]
pdf("Reactome_LvsNL.pdf", width=9, height=3.5)
dotplot(RAT05,showCategory=10)
dev.off()

##GSEA
library(msigdbr)
library(enrichplot)
gselist <- resLFCdata$log2FoldChange
names(gselist) <- resLFCdata$gene_name
gselist <- sort(gselist, decreasing = T)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, gene_symbol)
m_t2g$gs_name<- gsub("GO","",m_t2g$gs_name)
gsea <- GSEA(gselist, TERM2GENE = m_t2g)
write.csv(gsea@result,file="gsea_c5.csv")

gsea_select <- gsea@result[grep("NEUTROPHIL",gsea_1$Description),]
gsea_select$Description <- gsub("_"," ", gsea_select$Description)

for (i in seq(1:nrow(gsea_select))){
  label<- paste0("gsea_",gsea_select$ID[i],".pdf")
  pdf(label, width=9,height=8)  
  p <- gseaplot(gsea, geneSetID = gsea_select$ID[i], title = gsea_select$Description[i], ES_geom = "line",rel_heights = c(1.5, 0.5, 1), base_size = 21, subplots = 1:3, pvalue_table = T)
  print(p)
  dev.off()
} 
