#### fig S2:----
# fig S2A: umap of hSAT adipocytes colored by sample
p1 <- DimPlot(adipo5,group.by = "orig.ident", label=F,  pt.size = 0.5)
jpeg(paste0(output_path, "adipo5_bySample_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig S2B: umap of hVAT adipocytes colored by sample
p1 <- DimPlot(adipo10, group.by = "orig.ident", label=F,  pt.size = 0.5)
jpeg(paste0(output_path, "adipo10_bySample_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig S2C: umap of whole hSAT colored by sample
p1 <- DimPlot(sc5, group.by = "orig.ident", label=F,  pt.size = 0.5)
jpeg(paste0(output_path, "sc5_bySample_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig S2D: umap of whole hVAT colored by sample
p1 <- DimPlot(vis10, group.by = "orig.ident", label=F,  pt.size = 0.5)
jpeg(paste0(output_path, "vis10_bySample_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()
