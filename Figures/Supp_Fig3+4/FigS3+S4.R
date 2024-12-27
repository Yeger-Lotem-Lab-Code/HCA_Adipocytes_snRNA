#### fig S3:----
# big top20 DE for hSAT adipocytes heatmap
adipo5_logFC %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10_adipo5
top10_adipo5 <- top10_adipo5 %>% arrange(cluster)
p1 <- DoHeatmap(subset(adipo5, downsample=1000), features = top10_adipo5$gene, group.colors=cols_mapped_sc_adipo) +
  theme(text = element_text(size = 8))
jpeg(paste0(output_path, "adipo5_heatmap.jpeg"), width = 3000, height = 3400, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

#### fig S4:----
# big top20 DE for hVAT adipocytes heatmap
adipo10_logFC %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10_adipo10
top10_adipo10 <- top10_adipo10 %>% arrange(cluster)
p1 <- DoHeatmap(subset(adipo10, downsample=1000), features = top10_adipo10$gene, group.colors=cols_mapped_vis_adipo) +
  theme(text = element_text(size = 8))
jpeg(paste0(output_path, "adipo10_heatmap.jpeg"), width = 3000, height = 3400, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()
