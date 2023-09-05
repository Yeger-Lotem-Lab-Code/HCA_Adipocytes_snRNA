library(Seurat)
library(cowplot)
library(harmony)
library(ggpubr)

input_path <- ".../doubletFinder_analysis/doublets_removed/"
output_path <- ".../integration/"

# first create Seurat object of all the samples merged
sample_list <- list.files(input_path, pattern = "vis", full.names = F, recursive = F) #change "vis" according to tissue type
vis.list <- list()
for (i in 2:7) {
  print(sample_list[i])
  seurat_obj <- readRDS(file = (paste(input_path, sample_list[i], sep="")))
  vis.list <- c(vis.list, seurat_obj)
  rm(seurat_obj)
}

merg_vis <- merge(vis.list[[1]], y = vis.list[2:length(vis.list)], add.cell.ids = names(sample_list), project = "merg_vis")
merg_vis <- merg_vis %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(features = rownames(merg_vis), verbose = FALSE) %>% 
  RunPCA(verbose = FALSE)

options(repr.plot.height = 2.5, repr.plot.width = 6)
vis_harmony <- merg_vis %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
rm(merg_vis, vis.list, sample_list)

harmony_embeddings <- Embeddings(vis_harmony, 'harmony')
write.csv(harmony_embeddings, file = paste(output_path, "harmony_embeddings.csv", sep = ""))

p1 <- DimPlot(object = vis_harmony, reduction = "harmony", pt.size = 0.7, group.by = "orig.ident", do.return = TRUE)
p2 <- VlnPlot(object = vis_harmony, features = "harmony_1", group.by = "orig.ident", do.return = TRUE, pt.size = 0)
jpeg(paste(output_path, "harmony_PCs.jpeg", sep = ""), width = 1500, height = 1500)
plot_grid(p1,p2)
dev.off()
rm(p1, p2) 

vis_harmony <- vis_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 0.8)

p1 <- DimPlot(vis_harmony, reduction = "umap", group.by = "orig.ident", pt.size = 0.7)
p2 <- DimPlot(vis_harmony, reduction = "umap", label = TRUE, pt.size = 0.7)
jpeg(paste(output_path, "vis_harmony_bySample.jpeg", sep = ""), width = 1500, height = 1500)
print(p1)
dev.off()
jpeg(paste(output_path, "vis_harmony_byCluster.jpeg", sep = ""), width = 1500, height = 1500)
print(p2)
dev.off()
rm(p1,p2) 

saveRDS(vis_harmony, file = (paste(output_path, "vis_harmony.rds", sep="")))
