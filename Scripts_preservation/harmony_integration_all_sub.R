library(Seurat)
library(cowplot)
library(harmony)
library(ggpubr)
for (i in 7:10) {
  for (j in 1:10) {
    output_path <- paste0("/gpfs0/estiyl/users/tomul/sc_5/",i,"/",i,".",j,"/")
    
    merg11vis <- readRDS(file = paste0("/storage16/projects/Tom_Noa/sc_5/", i, "/", i, ".", j, "/sc_", i, ".", j, ".rds"))
    
    merg11vis <- merg11vis %>%
      Seurat::NormalizeData(verbose = FALSE) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
      ScaleData(features = rownames(merg11vis), verbose = FALSE) %>% 
      RunPCA(verbose = FALSE)
    
    options(repr.plot.height = 2.5, repr.plot.width = 6)
    vis11harmony <- merg11vis %>% 
      RunHarmony("orig.ident", plot_convergence = TRUE)
    # rm(merg11vis, vis.list, sample_list)
    rm(merg11vis)
    
    harmony_embeddings <- Embeddings(vis11harmony, 'harmony')
    write.csv(harmony_embeddings, file = paste(output_path, "harmony_embeddings.csv", sep = ""))
    
    p1 <- DimPlot(object = vis11harmony, reduction = "harmony", pt.size = 0.7, group.by = "orig.ident")
    p2 <- VlnPlot(object = vis11harmony, features = "harmony_1", group.by = "orig.ident", pt.size = 0)
    jpeg(paste(output_path, "harmony_PCs.jpeg", sep = ""), width = 1500, height = 1500)
    plot_grid(p1,p2)
    dev.off()
    rm(p1, p2) 
    
    vis11harmony <- vis11harmony %>% 
      RunUMAP(reduction = "harmony", dims = 1:50) %>% 
      FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
      FindClusters(resolution = 0.8)
    
    # changed the plot to umap_orig to compare it to the original umap dimentions
    p1 <- DimPlot(vis11harmony, reduction = "umap_orig", group.by = "orig.ident", pt.size = 0.7)
    p2 <- DimPlot(vis11harmony, reduction = "umap_orig", label = TRUE, pt.size = 0.7)
    jpeg(paste(output_path, "harmony_bySample_",i,".",j,".jpeg", sep = ""), width = 1500, height = 1500)
    print(p1)
    dev.off()
    jpeg(paste(output_path, "harmony_byCluster_",i,".",j,".jpeg", sep = ""), width = 1500, height = 1500)
    print(p2)
    dev.off()
    rm(p1,p2) 
    
    saveRDS(vis11harmony, file = (paste(output_path, "harmony_",i,".",j,".rds", sep="")))
    
  }
}

