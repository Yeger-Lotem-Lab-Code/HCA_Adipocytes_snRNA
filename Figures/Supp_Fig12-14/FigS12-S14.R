#### fig S12:----
# figure S12A: ForceAtlas2 hSAT ASPC+adipocytes colored by subtype, and their DPT pseudotime
p1 <- DimPlot(aspc_adipo5, reduction = "draw_graph_fa", group.by = "final_detailed" ,label=F,  pt.size = 0.7, cols = sc_aspc_adipo_col2)
jpeg(paste0(output_path, "aspc_adipo5_allcolored.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

#add pseudotime- from the code in Fig6

# figure S12B: umap of hSAT+hVAT ASPC integration colored by subtype and colored by tissue
p1 <- DimPlot(aspc_sc_vis, reduction = "umap", group.by = "final_detailed" ,label=F,  pt.size = 0.7, cols = aspc_sc_vis)
jpeg(paste0(output_path, "aspc_adipo5_allcolored.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- DimPlot(aspc_sc_vis, reduction = "umap", group.by = "Tissue", cols = c("#41ab5d","#dd3497"),label=F,  pt.size = 0.7)
jpeg(paste0(output_path, "aspc_adipo5_allcolored.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# figure S12C: ForceAtlas2 hVAT ASPC+adipocytes colored by subtype, and their DPT pseudotime
p1 <- DimPlot(aspc_adipo10, reduction = "draw_graph_fa", group.by = "final_detailed" ,label=F,  pt.size = 0.7, cols = vis_aspc_adipo_col2)
jpeg(paste0(output_path, "aspc_adipo10_allcolored.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

#add pseudotime-from the code in Fig6

#### fig S13:----
# figure S13A: umap hSAT ASPC+adipocytes colored by subtype after double filter
p1 <- DimPlot(aspc_adipo5_SD, reduction = "umap", group.by = "final_detailed" ,label=F,  pt.size = 0.7, cols = sc_aspc_adipo_col2)
jpeg(paste0(output_path, "aspc_adipo5_SD_UMAP_allcolored.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# figure S13B: violin plot of hSAT ASPC+adipocytes subtypes for "ADIPOQ", "AQP7", "LPL", "PLIN1"
features2 <- c("ADIPOQ","AQP7","LPL","PLIN1")
p1 <- VlnPlot(aspc_adipo5_SD, features=features2, stack = TRUE, flip = TRUE) +  theme(legend.position = "none")+
  scale_fill_manual(values=adipo_genes_colors)+
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)
jpeg(paste0(output_path, "aspc_adipo5_SD_violinPlot.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# figure S13C: ForceAtlas2 hSAT ASPC+adipocytes colored by subtype
p1 <- DimPlot(aspc_adipo5_SD, reduction = "draw_graph_fa", group.by = "final_detailed" ,label=F,  pt.size = 0.7, cols = sc_aspc_adipo_col2)
jpeg(paste0(output_path, "aspc_adipo5_SD_allcolored.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# figure S13D: umap hSAT ASPC+adipocytes colored by DPT pseudotime and palantir with different starting points
#add pseudotimes-from the code in Fig6

#### fig S14:----
# figure S14A: umap hVAT ASPC+adipocytes colored by subtype after double filter
p1 <- DimPlot(aspc_adipo10_SD, reduction = "umap", group.by = "final_detailed" ,label=F,  pt.size = 0.7, cols = vis_aspc_adipo_col2)
jpeg(paste0(output_path, "aspc_adipo10_SD_UMAP_allcolored.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# figure S14B: violin plot of hVAT ASPC+adipocytes subtypes for "ADIPOQ", "AQP7", "LPL", "PLIN1"
features2 <- c("ADIPOQ","AQP7","LPL","PLIN1")
p1 <- VlnPlot(aspc_adipo10_SD, features=features2, stack = TRUE, flip = TRUE) +  theme(legend.position = "none")+
  scale_fill_manual(values=adipo_genes_colors)+
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)
jpeg(paste0(output_path, "aspc_adipo10_SD_violinPlot.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# figure S14C: ForceAtlas2 hVAT ASPC+adipocytes colored by subtype
p1 <- DimPlot(aspc_adipo10_SD, reduction = "draw_graph_fa", group.by = "final_detailed" ,label=F,  pt.size = 0.7, cols = vis_aspc_adipo_col2)
jpeg(paste0(output_path, "aspc_adipo10_SD_allcolored.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# figure S14D: umap hVAT ASPC+adipocytes colored by DPT pseudotime and palantir with different starting points
#add pseudotimes-from the code in Fig6
