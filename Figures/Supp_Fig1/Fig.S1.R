#### fig S1: boxplots of paired samples:----
# boxplot for hSAT-hVAT adipocytes:
p1 <- paired_samples[paired_samples$cluster=='Adipocytes',] %>% ggplot(paired_samples, mapping=aes(x=Tissue, y=value, fill=Tissue))+ #fill=sex or fBMI
  geom_boxplot() + geom_boxplot(aes(color=Tissue),fatten = NULL, fill = NA) + theme_classic()+ scale_y_continuous(limits=c(0,70)) +
  scale_color_manual(values=c("#d95f0e","#fec44f"))+ scale_fill_manual(values=c("#d95f0e","#fec44f"))+
  ylab("Adipocytes %")+ theme_classic() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size = 18)) + 
  scale_x_discrete(labels=c("hSAT", "hVAT","Decon.\nhSAT", "Decon.\nhVAT"))

jpeg(paste0(output_path, "paired_adipocytes_boxplot.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# adipocytes: male and female
p1 <- paired_samples[paired_samples$cluster=='Adipocytes',] %>% ggplot(paired_samples, mapping=aes(x=Tissue, y=value, fill=Sex))+ #fill=sex or fBMI
  geom_boxplot() + geom_boxplot(aes(color=Sex),fatten = NULL, fill = NA, width=0.75) + theme_classic()+ scale_y_continuous(limits=c(0,70)) +
  scale_color_manual(values=sex_colors)+ scale_fill_manual(values=sex_colors)+
  ylab("Adipocytes %")+ theme_classic() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size = 18)) + 
  scale_x_discrete(labels=c("hSAT", "hVAT","Decon.\nhSAT", "Decon.\nhVAT"))

jpeg(paste0(output_path, "paired_adipocytes_boxplot_F-M.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# adipocytes; BMI group
p1 <- paired_samples[paired_samples$cluster=='Adipocytes',] %>% ggplot(paired_samples, mapping=aes(x=Tissue, y=value, fill=fBMI))+ #fill=sex or fBMI
  geom_boxplot(position = position_dodge(preserve = 'single')) + geom_boxplot(aes(color=fBMI),fatten = NULL, fill = NA, position = position_dodge(preserve = 'single')) + 
  theme_classic()+ scale_y_continuous(limits=c(0,70)) +
  scale_color_manual(values=BMI_colors)+ scale_fill_manual(values=BMI_colors)+
  ylab("Adipocytes %")+ theme_classic() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size = 18)) + 
  scale_x_discrete(labels=c("hSAT", "hVAT","Decon.\nhSAT", "Decon.\nhVAT"))

jpeg(paste0(output_path, "paired_adipocytes_boxplot_fBMI.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# adipocytes: age groups
p1 <- paired_samples[paired_samples$cluster=='Adipocytes',] %>% ggplot(paired_samples, mapping=aes(x=Tissue, y=value, fill=fage))+ #fill=sex or fBMI
  geom_boxplot(position = position_dodge(preserve = 'single')) + geom_boxplot(aes(color=fage),fatten = NULL, fill = NA,position = position_dodge(preserve = 'single')) + 
  theme_classic()+ scale_y_continuous(limits=c(0,70)) +
  scale_color_manual(values=age_colors)+ scale_fill_manual(values=age_colors)+
  ylab("Adipocytes %")+ theme_classic() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size = 18)) + 
  scale_x_discrete(labels=c("hSAT", "hVAT","Decon.\nhSAT", "Decon.\nhVAT"))

jpeg(paste0(output_path, "paired_adipocytes_boxplot_fage.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()
