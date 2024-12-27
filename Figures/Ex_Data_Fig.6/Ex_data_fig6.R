#### EXTENDED FIGURE 6:----
# ED FIG 6A: violin plot of hSAT+hVAT adipocytes for "ADIPOR1","ADIPOR2","CDH13","LEPR"
features1 <- c("ADIPOQ", "LEP","ADIPOR1","ADIPOR2","CDH13","LEPR")
p1 <- VlnPlot(adipo_vissc, group.by = "adipo_ann_new", features=features1, stack = TRUE, flip = TRUE) +  theme(legend.position = "none")+
  stat_summary(fun.y = median, geom='point', size = 6, colour = "black", shape = 95)
jpeg(paste0(output_path, "adipo_hSAT-hVAT_ADIPOQ_LEP_violinPlot.jpeg"), width = 1350, height = 1700, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# ED FIG 6B: violin plot of non adipocytes in hSAT+hVAT for "ADIPOR1","ADIPOR2","CDH13","LEPR"
features1 <- c("ADIPOR1","ADIPOR2","CDH13","LEPR")
p1 <- VlnPlot(sc5, group.by = "main_groups", features=features1, stack = TRUE, flip = TRUE) +  theme(legend.position = "none")+
  stat_summary(fun.y = median, geom='point', size = 6, colour = "black", shape = 95)
jpeg(paste0(output_path, "sc5_ADIPOQ_LEP_violinPlot.jpeg"), width = 1350, height = 1700, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- VlnPlot(vis10, group.by = "main_groups", features=features1, stack = TRUE, flip = TRUE) +  theme(legend.position = "none")+
  stat_summary(fun.y = median, geom='point', size = 6, colour = "black", shape = 95)
jpeg(paste0(output_path, "vis10_ADIPOQ_LEP_violinPlot.jpeg"), width = 1350, height = 1700, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()