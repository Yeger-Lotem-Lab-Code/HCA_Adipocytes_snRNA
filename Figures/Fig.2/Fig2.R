#### FIGURE 2:----
output_path <- "/final_figs/fig2/"

# fig 2A: umaps of sub-adipocytes+pieplots inset:
p1 <- DimPlot(adipo10, cols = cols_mapped_vis_adipo, label=F,  pt.size = 0.7)
jpeg(paste0(output_path, "adipo10_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- DimPlot(adipo5, cols = cols_mapped_sc_adipo, label=F,  pt.size = 0.7)
jpeg(paste0(output_path, "adipo5_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- ggplot(adipo10_disp, aes(x="", y=sum, fill=subtypes)) +
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=cols_mapped_vis_adipo)+ #
  theme(legend.position = "right", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16))+ 
  labs(fill = "Adipocytes\nsubtypes")
jpeg(paste0(output_path, "vis_adipo_subtypes_pieplot.jpeg"), width = 1900, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- ggplot(adipo5_disp, aes(x="", y=sum, fill=subtypes)) +
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=cols_mapped_sc_adipo)+ #
  theme(legend.position = "right", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16))+ 
  labs(fill = "Adipocytes\nsubtypes")
jpeg(paste0(output_path, "sc_adipo_subtypes_pieplot.jpeg"), width = 1900, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 2B: violin of adipocytes marker genes:
features2 <- c("ADIPOQ","AQP7","LPL","PLIN1")

p1 <- VlnPlot(adipo10, features=features2, stack = TRUE, flip = TRUE) +  theme(legend.position = "none")+
  scale_fill_manual(values=adipo_genes_colors)+
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)
jpeg(paste0(output_path, "adipo10_violinPlot.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- VlnPlot(adipo5, features=features2, stack = TRUE, flip = TRUE) +  theme(legend.position = "none")+
  scale_fill_manual(values=adipo_genes_colors)+
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)
jpeg(paste0(output_path, "adipo5_violinPlot.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)

# fig 2C: top20 marker genes for adipocytes: (too big, run on cluster)
adipo10_logFC %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10_adipo10
top10_adipo10 <- top10_adipo10 %>% arrange(cluster)
p1 <- DoHeatmap(subset(adipo10, downsample=1000), features = top10_adipo10$gene, group.colors=cols_mapped_vis_adipo) +
  theme(text = element_text(size = 8))
jpeg(paste0(output_path, "adipo10_heatmap.jpeg"), width = 3000, height = 3400, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

adipo5_logFC %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10_adipo5
top10_adipo5 <- top10_adipo5 %>% arrange(cluster)
p1 <- DoHeatmap(subset(adipo5, downsample=1000), features = top10_adipo5$gene, group.colors=cols_mapped_sc_adipo) +
  theme(text = element_text(size = 8))
jpeg(paste0(output_path, "adipo5_heatmap.jpeg"), width = 3000, height = 3400, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 2D: doubletFinder scores + MT-gene content:
v1 <- VlnPlot(adipo10, features = "percent.mt",pt.size = 0)+ ylim(0,20) + scale_fill_manual(values = cols_mapped_vis_adipo)+
  stat_summary(fun.y = median, geom='point', size = 6, colour = "black", shape = 95)
jpeg(paste0(output_path, "adipo10_percent.mt_vln.jpeg"), width = 3000, height = 3000, bg = "transparent",res=300)
print(v1)
dev.off()
v1 <- VlnPlot(adipo10, features = "DF.score.pANN",pt.size = 0)+ ylim(0,20) + scale_fill_manual(values = cols_mapped_vis_adipo)+
  stat_summary(fun.y = median, geom='point', size = 6, colour = "black", shape = 95)
jpeg(paste0(output_path, "adipo10_pDFpANN.mt_vln.jpeg"), width = 3000, height = 3000, bg = "transparent",res=300)
print(v1)
dev.off()

v1 <- VlnPlot(adipo5, features = "percent.mt",pt.size = 0)+ ylim(0,20) + scale_fill_manual(values = cols_mapped_sc_adipo)+
  stat_summary(fun.y = median, geom='point', size = 6, colour = "black", shape = 95)
jpeg(paste0(output_path, "adipo10_percent.mt_vln.jpeg"), width = 3000, height = 3000, bg = "transparent",res=300)
print(v1)
dev.off()
v1 <- VlnPlot(adipo5, features = "DF.score.pANN",pt.size = 0)+ ylim(0,20) + scale_fill_manual(values = cols_mapped_sc_adipo)+
  stat_summary(fun.y = median, geom='point', size = 6, colour = "black", shape = 95)
jpeg(paste0(output_path, "adipo10_pDFpANN.mt_vln.jpeg"), width = 3000, height = 3000, bg = "transparent",res=300)
print(v1)
dev.off()

# fig 2E: UMAP of hSAT+hVAT integrated:
p1 <- DimPlot(adipo_vissc, group.by = "Tissue", cols = c("#41ab5d","#dd3497"),  label=F,  pt.size = 0.2)
jpeg(paste0(output_path, "adipo_vissc_tiss_UMAP_noLabel.jpeg"), width = 2100, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# fig 2F: barplot for hSAT-hVAT in integration:
p1 <- ggplot(adipo_vissc@meta.data, aes(x=com_idents, fill=adipo_ann_new)) + geom_bar(position = "fill")+ 
  scale_fill_manual(values=cols4bars)+
  ggtitle("SAT and VAT adipocytes proportions") +
  xlab("Integrated subtypes") + ylab("Proportions") + theme_classic() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, size=12),
        plot.title = element_text(hjust = 0.5, size = 14)) + 
  labs(fill = "Adipocytes\nsubtypes")+ 
  scale_x_discrete(labels=c("iASV1\n22,571", "iASV2\n431", "iASV3\n837","iASV4\n440", 
                            "iASV5\n3,064", "iASV6\n222", "iASV7\n100")) 
jpeg(paste0(output_path, "SAT-VAT_integrated_subtypes_barplot.jpeg"), width = 1700, height = 1500, quality = 100, res=300)
print(p1)
dev.off()

# fig 2G: table: taken from table S7A
