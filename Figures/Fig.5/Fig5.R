
#### FIGURE 5:----
output_path <- "/final_figs/fig5/"

# fig 5A: umap of adipo5+adipo10 for PTPRB:
p1 <- FeaturePlot(adipo5, features = "PTPRB", label=F,  pt.size = 0.6,order = T)
jpeg(paste0(output_path, "adipo5_PTPRB.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(adipo10, features = "PTPRB", label=F,  pt.size = 0.6,order = T)
jpeg(paste0(output_path, "adipo10_PTPRB.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 5B: umap of adipo5+adipo10 for PDE4D:
p1 <- FeaturePlot(adipo5, features = "PDE4D", label=F,  pt.size = 0.6,order = T)
jpeg(paste0(output_path, "adipo5_PDE4D.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(adipo10, features = "PDE4D", label=F,  pt.size = 0.6,order = T)
jpeg(paste0(output_path, "adipo10_PDE4D.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 5C: umap of adipo5+adipo10 for SKAP1:
p1 <- FeaturePlot(adipo5, features = "SKAP1", label=F,  pt.size = 0.6,order = T)
jpeg(paste0(output_path, "adipo5_SKAP1.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(adipo10, features = "SKAP1", label=F, pt.size = 0.6,order = T)
jpeg(paste0(output_path, "adipo10_SKAP1.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 5D: umap of adipo5+adipo10 for ANK2:
p1 <- FeaturePlot(adipo5, features = "ANK2", label=F,  pt.size = 0.6,order = T)
jpeg(paste0(output_path, "adipo5_ANK2.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(adipo10, features = "ANK2", label=F,  pt.size = 0.6,order = T)
jpeg(paste0(output_path, "adipo10_ANK2.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 5I: ridgeplot for nuc-MT% in hVAT adipocytes:
p1 <-  ggplot(adipo10_meta, aes(x = percent.nuc_mito , y = names,  fill=names)) +  
  scale_y_discrete(limits = rev) +
  geom_density_ridges(scale = 0.9, alpha=0.7) +
  theme_ridges(center_axis_labels = TRUE) +
  scale_x_continuous(limits = c(0, 10)) +  # Set x-axis range from 0 to 5
  scale_fill_manual(values = cols_mapped_vis_adipo) 

jpeg(paste0(output_path, "adipo10_nucMT_ridge.jpeg"), width = 3000, height = 3000, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 5J: umap of adipo10 for MT-CO2:
p1 <- FeaturePlot(adipo10, features = "MT-CO2", label=F,  pt.size = 0.6,order = T)
jpeg(paste0(output_path, "adipo10_MT-CO2.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()


