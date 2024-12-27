#### fig S9:----
# figure S9A:
p1 <- ggplot(ML_results, mapping=aes(x=var_tis, y=Val, fill=Var))+ 
  geom_boxplot() + geom_boxplot(aes(color=Var),fatten = NULL, fill = NA,position = position_dodge(preserve = 'single')) + theme_classic()+ scale_y_continuous(limits=c(0.65,1)) +
  scale_color_manual(values=c("#7570b3","#d95f02","#1b9e77"))+
  scale_fill_manual(values=c("#7570b3","#d95f02","#1b9e77"))+
  ylab("")+xlab("")+
  geom_jitter(color="black", size=1, alpha=0.9, width = 0.2) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, size=10),
        axis.text.y = element_text(size=18))

jpeg(paste0(output_path, "AUC_PR_ML_results_boxplot.jpeg"), width = 1800, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# figure S9B:
# vis:
vis_sum <- as.data.frame(table(vis_predict$predicted))
p1 <- ggplot(vis_sum, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=cols_mapped_vis_adipo)+ #,labels = c("Specialized adipocytes", "Classical adipocytes")
  theme(legend.position = "right", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16))+ 
  labs(fill = "Adipocytes\nsubtypes")
jpeg(paste0(output_path, "vis_predicted_pieplot_0.2.jpeg"), width = 1900, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# sc:
sc_sum <- as.data.frame(table(sc_predict$predicted))
p1 <- ggplot(sc_sum, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=cols_mapped_sc_adipo)+ #,labels = c("Specialized adipocytes", "Classical adipocytes")
  theme(legend.position = "right", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16))+ 
  labs(fill = "Adipocytes\nsubtypes")
jpeg(paste0(output_path, "sub_predicted_pieplot_0.2.jpeg"), width = 1900, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# figure S9C: boxplot of adipocytes as 100%, classical+non-classical in each tissue
p1 <- ggplot(cVs_adipo, mapping=aes(x=Tissue, y=value_new, fill=cluster))+ 
  geom_boxplot() + 
  geom_boxplot(aes(color=cluster),fatten = NULL, fill = NA) + 
  theme_classic()+ scale_y_continuous(limits=c(0,100)) +
  scale_color_manual(values=c(cVs_colors,cVs_colors))+
  scale_fill_manual(values=c(cVs_colors,cVs_colors))+ 
  xlab("Tissue") + 
  ylab("Adipocytes %")+ theme_classic() +
  geom_jitter(color="black", size=0.7, alpha=0.9) +
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size = 18)) + 
  scale_x_discrete(labels=c("hSAT", "hVAT","Decon.\nhSAT", "Decon.\nhVAT"))

jpeg(paste0(output_path, "cVs_as_100_decon_step1_tissue_boxplot.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()
