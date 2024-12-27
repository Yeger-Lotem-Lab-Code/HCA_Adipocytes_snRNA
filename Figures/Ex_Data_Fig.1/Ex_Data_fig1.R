#### EXTENDED FIGURE 1:----
output_path <- "/final_figs/extfig/"
# 1 big umap of whols hVAT colored by cell type, 6 small umaps of 
# "PTPRC" , "PDE4D" , "ANK2" , "PTPRB" , "SKAP1" , "MT-CO2"
p1 <- DimPlot(vis10, cols = cols_mapped, label=F,  pt.size = 0.5)
jpeg(paste0(output_path, "vis10_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(vis10, features = "PTPRC", label=F,  pt.size = 0.5)
jpeg(paste0(output_path, "vis10_PTPRC_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(vis10, features = "PDE4D", label=F,  pt.size = 0.5)
jpeg(paste0(output_path, "vis10_PDE4D_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(vis10, features = "ANK2", label=F,  pt.size = 0.5)
jpeg(paste0(output_path, "vis10_ANK2_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(vis10, features = "PTPRB", label=F,  pt.size = 0.5)
jpeg(paste0(output_path, "vis10_PTPRB_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(vis10, features = "SKAP1", label=F,  pt.size = 0.5)
jpeg(paste0(output_path, "vis10_SKAP1_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(vis10, features = "MT-CO2", label=F,  pt.size = 0.5)
jpeg(paste0(output_path, "vis10_MT-CO2_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()