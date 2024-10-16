library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGREST)

input_path <- ".../DE/"
output_path <- ".../DE/GSEA/"

# for gseGO/KEGG
FC_table <- read.csv(paste0(input_path,"logFC_markers.csv"),row.names = 1, header = T)
FC_table <- FC_table[FC_table$pct.1>=0.1,]
FC_table <- FC_table[order(FC_table$avg_log2FC, decreasing = TRUE),]  

FC_table$cluster <- factor(FC_table$cluster)
clusters <- levels(FC_table$cluster)

for (i in 1:length(clusters)) {
  geneList <- FC_table[FC_table$cluster==i, c("gene", "avg_log2FC")]
  geneList <- geneList[order(geneList$avg_log2FC, decreasing = TRUE),] 
  ranked_geneList <- geneList$avg_log2FC
  names(ranked_geneList) <- geneList$gene
  
  symbol2entrezid <- bitr(geneList$gene, fromType="SYMBOL", toType="ENTREZID", 
                          OrgDb="org.Hs.eg.db")
  symbol2entrezid <- symbol2entrezid %>% distinct(SYMBOL,.keep_all = TRUE)
  
  geneList_KEGG <- geneList %>% filter(gene %in% symbol2entrezid$SYMBOL)
  ranked4KEGG <- geneList_KEGG$avg_log2FC
  names(ranked4KEGG) <- symbol2entrezid$ENTREZID

  GO_results <- gseGO(geneList = ranked_geneList, keyType = "SYMBOL",
                      OrgDb = org.Hs.eg.db, ont = "ALL", minGSSize = 2)
  write.csv(GO_results, file = paste0(output_path, "cluster_", i, "_gseGO.csv"))
  
  KEGG_results <- gseKEGG(geneList = ranked4KEGG, minGSSize = 2, organism = "hsa") # has to use ENTREZID, no symbol
  KEGG_results <-  setReadable(KEGG_results, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  write.csv(KEGG_results, file = paste0(output_path, "cluster_", i, "_gseKEGG.csv"))
  
  rm(geneList, ranked_geneList, symbol2entrezid, GO_results,ranked4KEGG, geneList_KEGG, KEGG_results)
}

