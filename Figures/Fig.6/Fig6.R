#### FIGURE 6:----
output_path <- "/final_figs/fig6/"

# fig 6A: umap of hSAT PDGFRA
p1 <- FeaturePlot(sc5, features = "PDGFRA" ,label=F,  pt.size = 0.2)
jpeg(paste0(output_path, "aspc_adipo5_PDGFRA.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 6B: pieplot for hSAT ASPCs:
p1 <- ggplot(sc_aspc_disp, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=aspc_only)+ #
  theme(legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))+
  labs(fill = "ASPCs\nsubtypes")
jpeg(paste0(output_path, "sc_aspc_subtypes_pieplot.jpeg"), width = 1900, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 6C: vlnplot for hSAT ASPCs marker:
features4ASPC <- c("PDGFRA", "DPP4", "CD55", "PPARG", "ITGB1", "THY1", "CD9", "PDE4D","ALDH1A3", "WT1","F3")
p1 <- VlnPlot(aspc_sc_only, features=features4ASPC, stack = TRUE, flip = TRUE) +  theme(legend.position = "none")+
  scale_y_continuous(position="left")+
  stat_summary(fun.y = median, geom='point', size = 6, colour = "black", shape = 95)
jpeg(paste0(output_path, "aspc_adipo5_markers_violinPlot.jpeg"), width = 1350, height = 1700, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 6D: force-map for hSAT ASPC->non-classical adipo->classical adipo (with inset of ASPC subtypes):
# force-map colored by cell types:
p1 <- DimPlot(aspc_adipo5, reduction = "draw_graph_fa", group.by = "final_detailed" ,label=F,  pt.size = 0.7, cols = sc_aspc_adipo_col)
jpeg(paste0(output_path, "aspc_adipo5.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()
# only ASPC colored:
p1 <- DimPlot(aspc_sc_only, reduction = "draw_graph_fa", group.by = "final_detailed" ,label=F,  pt.size = 0.7, cols = aspc_only)
jpeg(paste0(output_path, "aspc_sc_only.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()
# force-map colored by pseudotime score:
#dim reduction= draw_graph_fa

# fig 6E: umap of hVAT PDGFRA
p1 <- FeaturePlot(vis10, features = "PDGFRA", label=F,  pt.size = 0.2)
jpeg(paste0(output_path, "aspc_adipo10_PDGFRA.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 6F: pieplot for hVAT ASPCs:
p1 <- ggplot(vis_aspc_disp, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=aspc_only)+ #
  theme(legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))+
  labs(fill = "ASPCs\nsubtypes")
jpeg(paste0(output_path, "vis_aspc_subtypes_pieplot.jpeg"), width = 1900, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 6G: vlnplot for hVAT ASPCs marker:
features4ASPC <- c("PDGFRA", "DPP4", "CD55", "PPARG", "ITGB1", "THY1", "CD9", "PDE4D","ALDH1A3", "WT1","F3")
p1 <- VlnPlot(aspc_vis_only, features=features4ASPC, stack = TRUE, flip = TRUE) +  theme(legend.position = "none")+
  scale_y_continuous(position="left")+
  stat_summary(fun.y = median, geom='point', size = 6, colour = "black", shape = 95)
jpeg(paste0(output_path, "aspc_adipo10_markers_violinPlot.jpeg"), width = 1350, height = 1700, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 6H: force-map for hVAT ASPC->non-classical adipo->classical adipo (with inset of ASPC subtypes):
# force-map colored by cell types:
p1 <- DimPlot(aspc_adipo10, reduction = "draw_graph_fa", group.by = "final_detailed" ,label=F,  pt.size = 0.7, cols = vis_aspc_adipo_col)
jpeg(paste0(output_path, "aspc_adipo10.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()
# only ASPC colored:
p1 <- DimPlot(aspc_vis_only, reduction = "draw_graph_fa", group.by = "final_detailed" ,label=F,  pt.size = 0.7, cols = aspc_only)
jpeg(paste0(output_path, "aspc_vis_only.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()
# force-map colored by pseudotime score:
