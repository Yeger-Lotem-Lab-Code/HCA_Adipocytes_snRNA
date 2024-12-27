#### fig S5:----
# fig S5A: umap of double filtered hSAT adipocytes, colored by OG subtypes
p1 <- DimPlot(adipo5_scrub_DF, cols = cols_mapped_sc_adipo, label=F,  pt.size = 1)
jpeg(paste0(output_path, "adipo5_scrub_DF_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig S5B: umap of double filtered hVAT adipocytes, colored by OG subtypes
p1 <- DimPlot(adipo10_scrub_DF, cols = cols_mapped_vis_adipo, label=F,  pt.size = 1)
jpeg(paste0(output_path, "adipo10_scrub_DF_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

#### fig S6:----
# fig S6A: umap of Emont-pipeline hSAT adipocytes, colored by in-house subtypes
p1 <- DimPlot(adipo5_SCT, cols = cols_mapped_sc_adipo, label=F,  pt.size = 1)
jpeg(paste0(output_path, "adipo5_SCT_UMAP_SA_names.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig S6B: umap of Emont-pipeline hSAT adipocytes, colored by clusters
p1 <- DimPlot(adipo5_SCT,group.by="seurat_clusters", label=F,  pt.size = 1)
jpeg(paste0(output_path, "adipo5_SCT_UMAP_clusters.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig S6C: umap of Emont-pipeline hVAT adipocytes, colored by in-house subtypes
p1 <- DimPlot(adipo10_SCT, cols = cols_mapped_vis_adipo, label=F,  pt.size = 1)
jpeg(paste0(output_path, "adipo10_SCT_UMAP_VA_names.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig S6D: umap of Emont-pipeline hVAT adipocytes, colored by clusters
p1 <- DimPlot(adipo10_SCT,group.by="seurat_clusters", label=F,  pt.size = 1)
jpeg(paste0(output_path, "adipo10_SCT_UMAP_clusters.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

#### fig S7:----
# fig S7A: umap of in-house+Emont hSAT adipocytes integration, colored by in-house subtypes
p1 <- DimPlot(sc_house_emont, cols = cols_mapped_2emont_sc_adipo, label=F,  pt.size = 1)
jpeg(paste0(output_path, "house_emont_hSAT_UMAP_bySAnames.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig S7B: umap of in-house+Emont hSAT adipocytes integration, colored by Emont subtypes
p1 <- DimPlot(sc_house_emont, cols = cols_mapped_1emont_sc_adipo, label=F,  pt.size = 1)
jpeg(paste0(output_path, "house_emont_hSAT_UMAP_byEmontnames.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig S7C: umap of in-house+Emont hVAT adipocytes integration, colored by in-house subtypes
p1 <- DimPlot(vis_house_emont, cols = cols_mapped_2emont_vis_adipo, label=F,  pt.size = 1)
jpeg(paste0(output_path, "house_emont_hVAT_UMAP_byVAnames.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig S7D: umap of in-house+Emont hVAT adipocytes integration, colored by Emont subtypes
p1 <- DimPlot(vis_house_emont, cols = cols_mapped_1emont_vis_adipo, label=F,  pt.size = 1)
jpeg(paste0(output_path, "house_emont_hVAT_UMAP_byEmontnames.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()