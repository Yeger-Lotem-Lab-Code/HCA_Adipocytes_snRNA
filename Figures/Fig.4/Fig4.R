# scripts for each main fig:
#### FIGURE 4:----
output_path <- "/final_figs/fig4/"

# fig 4A: scenic summary of selected regulon activity:
tab4fig <- read.csv(paste0(input_path, "scenic table for figure.csv")) # table of chosen regulons adaprwd from table S9
long <- tab4fig %>% 
  pivot_longer(
    cols = `SA1`:`VA8`, 
    names_to = "subset",
    values_to = "percentage"
  )
long$X <- factor(long$X, levels=c("ATF2","CEBPB","CEBPA","ATOH8","MECOM",
                                  "MEF2C","RELB","BACH1","RUNX1","RELA",
                                  "STAT4","STAT6","GLIS3","GLI2","PRRX1",
                                  "MYC","GABPB1","GABPA"))
plot <- ggplot(data = long, mapping = aes(x = subset, y = X, colour = percentage)) +
  geom_point(aes(size = after_stat(long$percentage), group = 1)) +
  scale_size_area(name="percentage",max_size = 10) + scale_y_discrete(limits=rev) +
  scale_colour_gradientn(name="Percentage",limits = c(0,100),
                         colours=c("#ffffcc", "#fd8d3c", "#800026"),
                         breaks=c(0,20,40,60,80,100)) +
  theme_classic()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=16),
        axis.text.y = element_text(angle = 0, size=18))
jpeg(paste0(input_path, "regulon_summary.jpeg"), width = 3000, height = 3000, quality = 100, bg = "transparent",res=300)
print(plot)
dev.off()

# fig 4B: umap integration of in-house and Emont et al. hSAT adipocytes + barplot:
p1 <- DimPlot(sc_house_emont, group.by = "data_source", cols = c("#969696","#41ab5d"),  label=F,  pt.size = 0.2)
jpeg(paste0(output_path, "sc_house_emont_tiss_UMAP_noLabel.jpeg"), width = 2100, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- ggplot(sc_house_emont@meta.data, aes(x=com_idents, fill=barplot)) + geom_bar(position = "fill")+ 
  scale_fill_manual(values=cols4bars_sc_emont2)+
  ggtitle("SAT and Emont adipocytes proportions") +
  xlab("Integrated subtypes") + ylab("Proportions") + theme_classic() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, size=12),
        plot.title = element_text(hjust = 0.5, size = 14)) + 
  labs(fill = "Adipocytes\nsubtypes")+ 
  scale_x_discrete(labels=c("i-SA1\n24,062", "i-SA2\n233", "i-SA3\n514","i-SA4\n279", "i-SA5\n1,296", "i-SA6\n139", "i-SA7\n78")) 
jpeg(paste0(output_path, "emont-SAT_integrated_subtypes_barplot.jpeg"), width = 1700, height = 1500, quality = 100, res=300)
print(p1)
dev.off()

# fig 4C: umap integration of in-house and Emont et al. hVAT adipocytes + barplot:
p1 <- DimPlot(vis_house_emont, group.by = "data_source", cols = c("#969696","#dd3497"),  label=F,  pt.size = 0.2)
jpeg(paste0(output_path, "vis_house_emont_tiss_UMAP_noLabel.jpeg"), width = 2100, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- ggplot(vis_house_emont@meta.data, aes(x=new_clusters, fill=barplot)) + geom_bar(position = "fill")+ 
  scale_fill_manual(values=cols4bars_vis_emont2)+
  ggtitle("VAT and Emont adipocytes proportions") +
  xlab("Integrated subtypes") + ylab("Proportions") + theme_classic() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, size=12),
        plot.title = element_text(hjust = 0.5, size = 14)) + 
  labs(fill = "Adipocytes\nsubtypes")+ 
  scale_x_discrete(labels=c("i-VA1\n22,834", "i-VA2\n346", "i-VA3\n559","i-VA4\n277", "i-VA5\n2,919")) 
jpeg(paste0(output_path, "emont-VAT_integrated_subtypes_barplot.jpeg"), width = 1700, height = 1500, quality = 100, res=300)
print(p1)
dev.off()

# pieplots of specialized subtypes in in-house hSAT + inetegrated in-house+Emont:
p1 <- ggplot(spec5, aes(x="", y=sum, fill=subtypes)) +
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=cols_mapped_sc_adipo)+ #
  theme(legend.position = "right", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16))+ 
  labs(fill = "Adipocytes\nsubtypes")
jpeg(paste0(output_path, "sc_adipo_special_only_pieplot.jpeg"), width = 1900, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

sc_he_spec <- sc_house_emont_disp[-1,]
p1 <- ggplot(sc_he_spec, aes(x="", y=emont_sum, fill=iVAT)) +
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=cols_mapped_iSAT)+ #
  theme(legend.position = "right", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16))+ 
  labs(fill = "Integrated SAT\nadipocytes")
jpeg(paste0(input_path, "i-SAT_special_pieplot.jpeg"), width = 1900, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# pieplots of specialized subtypes in in-house hVAT + inetegrated in-house+Emont:
p1 <- ggplot(spec10, aes(x="", y=sum, fill=subtypes)) +
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=cols_mapped_vis_adipo)+ #
  theme(legend.position = "right", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16))+ 
  labs(fill = "Adipocytes\nsubtypes")
jpeg(paste0(output_path, "vis_adipo_special_only_pieplot.jpeg"), width = 1900, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

vis_he_spec <- vis_house_emont_disp[-2,]
p1 <- ggplot(vis_he_spec, aes(x="", y=emont_sum, fill=iVAT)) +
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=cols_mapped_iVAT)+ #
  theme(legend.position = "right", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16))+ 
  labs(fill = "Integrated VAT\nadipocytes")
jpeg(paste0(output_path, "i-VAT_special_pieplot.jpeg"), width = 1900, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()


