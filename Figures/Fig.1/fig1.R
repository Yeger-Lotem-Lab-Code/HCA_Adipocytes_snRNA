# scripts for each main fig:

# libraries:

#### colors:----
cols_mapped <- c("Adipocytes" = "#ff7f00",
                 "ASPC" = "#fdbf6f" ,
                 "B cells" = "#c51b7d", #
                 "Dendritic cells" = "#a6cee3", #
                 "Endothelial" = "#6a3d9a" ,
                 "Endothelial (lymph.)" = "#cab2d6" ,
                 "Macrophages" = "#33a02c",
                 "Mast cells" = "#e5e555",
                 "Mesothelial" = "#e31a1c",
                 "Monocytes" =  "#1f78b4",
                 "Pericytes/SMC" = "#35978f" , 
                 "Pre-mesothelial" = "#fb9a99",#
                 "Pre-myeloids" = "#b2e485" ,
                 "T cells" = "#b15928",
                 "Unknown" = "#969696")

cols_mapped_vis_adipo <- c("VA1" = "#ff7f00", "VA2" = "#1f78b4", "VA3" = "#33a02c",
                           "VA4" = "#b2df8a", "VA5" = "#6a3d9a", "VA6" = "#a6cee3", 
                           "VA7" = "#cab2d6","VA8" = "#fdbf6f")

cols_mapped_iVAT <- c("iEVA1" = "#ff7f00", "iEVA2" = "#1f78b4", "iEVA3" = "#33a02c",
                      "iEVA4" = "#b2df8a", "iEVA5" = "#6a3d9a")

cols_mapped_sc_adipo <- c("SA1" = "#ff7f00", "SA2" = "#1f78b4", "SA3" = "#33a02c",
                          "SA4" = "#b2df8a", "SA5" = "#6a3d9a", "SA6" = "#fb9a99", "SA7" = "#e31a1c")

cols_mapped_iSAT <- c("iESA1" = "#ff7f00", "iESA2" = "#1f78b4", "iESA3" = "#33a02c",
                      "iESA4" = "#b2df8a", "iESA5" = "#6a3d9a", "iESA6" = "#fb9a99", "iESA7" = "#e31a1c")

cols_mapped_vissc <- c("Classical" = "#ff7f00", "Special 2" = "#1f78b4","Special 3" = "#33a02c", 
                       "Special 4" = "#b2df8a", "Special 5" = "#6a3d9a", "Special 6"="#fb9a99",
                       "Special V7" ="#cab2d6" , "Special S7V5" = "#e31a1c", "Special 7" ="#e31a1c" )

cols_mapped_2emont_vis_adipo <- c("VA1" = "#ff7f00", "VA2" = "#1f78b4", "VA3" = "#33a02c",
                                  "VA4" = "#b2df8a", "VA5" = "#6a3d9a", "VA6" = "#a6cee3", 
                                  "VA7" = "#cab2d6","VA8" = "#fdbf6f",
                                  "hAd1" = "#969696", "hAd2" = "#969696", "hAd3" = "#969696", 
                                  "hAd4" = "#969696", "hAd5" = "#969696", "hAd6" = "#969696", 
                                  "hAd7" = "#969696")

cols_mapped_1emont_vis_adipo <- c("hAd1" = "#ff7f00", "hAd2" = "#1f78b4", "hAd3" = "#33a02c",
                                  "hAd4" = "#b2df8a", "hAd5" = "#6a3d9a", "hAd6" = "#a6cee3", 
                                  "hAd7" = "#cab2d6",
                                  "VA1" = "#969696", "VA2" = "#969696", "VA3" = "#969696",
                                  "VA4" = "#969696", "VA5" = "#969696", "VA6" = "#969696", 
                                  "VA7" = "#969696","VA8" = "#969696")

cols_mapped_1emont_sc_adipo <- c("hAd1" = "#ff7f00", "hAd2" = "#1f78b4", "hAd3" = "#33a02c","hAd4" = "#b2df8a", 
                                 "hAd5" = "#6a3d9a", "hAd6" = "#fb9a99", "hAd7" = "#e31a1c",
                                 "SA1" = "#969696", "SA2" = "#969696", "SA3" = "#969696", 
                                 "SA4" = "#969696", "SA5" = "#969696", "SA6" = "#969696", 
                                 "SA7" = "#969696")

cols_mapped_2emont_sc_adipo <- c("SA1" = "#ff7f00", "SA2" = "#1f78b4", "SA3" = "#33a02c","SA4" = "#b2df8a", 
                                 "SA5" = "#6a3d9a", "SA6" = "#fb9a99", "SA7" = "#e31a1c",
                                 "hAd1" = "#969696", "hAd2" = "#969696", "hAd3" = "#969696", 
                                 "hAd4" = "#969696", "hAd5" = "#969696", "hAd6" = "#969696", 
                                 "hAd7" = "#969696")


tis_colors = c("#d95f0e","#fec44f")
tis_colors_decon = c("#d95f0e","#fec44f","#999999", "#CCCCCC")
cVs_colors <- c("#ff7f00","#6a3d9a","#fdbf6f","#cab2d6")
sex_colors = c("#b2df8a", "#33a02c")
cVs_sex_colors <- c("#b2df8a","#b2df8a", "#33a02c","#33a02c")
BMI_colors = c("#deebf7","#9ecae1","#3182bd")
age_colors = c("#cbc9e2","#9e9ac8","#756bb1","#54278f")
adipo_genes_colors = c("#fed976","#feb24c","#fd8d3c","#f03b20")

#### table preparations:----
# in-house meta-data:
input_path <- "C:/Users/or711/OneDrive - post.bgu.ac.il/MSc-PhD/adipo/meta_data/"
house_meta <-  read.csv(paste0(input_path, "clinical_parameters_2024-10-20.csv"), header=TRUE)
house_meta$Age <- round(house_meta$Age, digits = 1)
# 72 bulk samples meta-data:
input_path <- "C:/Users/or711/OneDrive - post.bgu.ac.il/MSc-PhD/adipo/deconvolution/"
meta72 <- readRDS(paste0(input_path, "Phenotypes_filtered_downsampled_additionalPhenos.rds"))
meta72 <- meta72 %>% mutate(fBMI = case_when(BMI <= 30 ~ "A", 
                                             BMI > 30 & BMI < 40 ~ "B",
                                             BMI >= 40 ~ "C"))
meta72 <- meta72 %>% mutate(fage = case_when(Age < 35 ~ "A", 
                                             Age >= 35 & Age < 50 ~ "B",
                                             Age >= 50 & Age < 65 ~ "C",
                                             Age >= 65 ~ "D" ))
meta_step1 <- meta72[meta72$Step=='1',]
meta_step2 <- meta72[meta72$Step=='2',]

# deconvolution prediction results: 2 .txt files for each model (CvS, main-names - only look at aggregated results (all endo types together))
# final decon+house cell%:
percentage_concat_CvS <- "agg_cVs_summarized data table.csv"
percentage_concat_comnames <- "agg_com_names_summarized data table.csv"


# seurat objects:----
vis10 <- readRDS(paste0(obj_path, "vis10_diet.rds"))
sc5 <- readRDS(paste0(obj_path, "sc5_diet.rds"))
adipo10 <- readRDS(paste0(obj_path, "adipo10_diet.rds"))
adipo5 <- readRDS(paste0(obj_path, "adipo5_diet.rds"))


#### FIGURE 1:----
output_path <- "C:/Users/or711/OneDrive - post.bgu.ac.il/MSc-PhD/adipo/figs4paper/final_figs/fig1/"

# figure 1A: pieplots of main clinical parameters:
# BMI: 3 colors (blues)
BMI_disp <- data.frame(with(house_meta, table(fBMI)))
p1 <- ggplot(BMI_disp, aes(x="", y=Freq, fill=fBMI)) +
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=BMI_colors, labels = c("BMI<=30", "30<BMI<40", "40<=BMI"))+
  theme(legend.position = "bottom", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16))+ 
  labs(fill = "BMI groups")
jpeg(paste0(output_path, "BMI_pieplot.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# sex: 2 (green)
sex_disp <- data.frame(with(house_meta, table(Sex)))
p1 <- ggplot(sex_disp, aes(x="", y=Freq, fill=Sex)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=sex_colors, labels = c("Female", "Male"))+#
  theme(legend.position = "bottom", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16))+ 
  labs(fill = "Sex")
jpeg(paste0(output_path, "decon_sex_pieplot.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# age: 4 (purple)
age_disp <- data.frame(with(house_meta, table(fage)))
p1 <- ggplot(age_disp, aes(x="", y=Freq, fill=fage)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=age_colors,labels = c("20<=Age<35", "35<=Age<50","50<=Age<65","65<=Age"))+#
  theme(legend.position = "bottom", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16))+ 
  labs(fill = "Age groups")
jpeg(paste0(output_path, "age_pieplot.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# tissues: 2 colors (orange)
tis_disp <- data.frame(Freq=c(5, 10), tissue=c("Subcutaneous", "Visceral"))
p1 <- ggplot(tis_disp, aes(x="", y=Freq, fill=tissue)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=tis_colors)+#,labels = c("", "35<=Age<50","50<=Age<65","65<=Age")
  theme(legend.position = "bottom", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16))+ 
  labs(fill = "Tissue type")
jpeg(paste0(output_path, "tissue_pieplot.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# figure 1B: UMAPS of adiponectin exp. without labels:
p1 <- FeaturePlot(vis10_diet, features = "ADIPOQ", label=F,  pt.size = 0.2)
jpeg(paste0(output_path, "vis10_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(sc5_diet, features = "ADIPOQ", label=F,  pt.size = 0.2)
jpeg(paste0(output_path, "sc5_UMAP_noLabel.jpeg"), width = 1750, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# figure 1C: violinplot of marker genes:
features1 <- c("ADIPOQ","AQP7","LPL","PLIN1","PDGFRA", "MS4A1","CLEC9A","VWF","FLT4","CD163","KIT","MSLN","ITGAM","TAGLN","DPP4","PROS1","CD247")
vis10_diet_sub <- subset(vis10_diet, idents = "Unknown", invert=T)
p1 <- VlnPlot(vis10_diet_sub, features=features1, stack = TRUE, flip = TRUE) +  theme(legend.position = "none")+
  scale_y_continuous(position="right")+
  stat_summary(fun.y = median, geom='point', size = 6, colour = "black", shape = 95)
jpeg(paste0(output_path, "vis_violinPlot.jpeg"), width = 1350, height = 1700, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

p2 <- VlnPlot(sc5_diet, features=features1, stack = TRUE, flip = TRUE) +  theme(legend.position = "none")+
  stat_summary(fun.y = median, geom='point', size = 6, colour = "black", shape = 95)
jpeg(paste0(output_path, "sc_violinPlot.jpeg"), width = 1000, height = 1700, quality = 100, bg = "transparent",res=300)
print(p2)
dev.off()

# figure 1D: boxplots of %adipocytes by clinical parameters;
samples <- c("3382", "3384", "3387", "3399", "3427")
paired_adipocytes <- percentage_concat_comnames[percentage_concat_comnames$variable %in% samples,]
# regular vis and sc:
p1 <- percentage_concat_comnames[percentage_concat_comnames$cluster=='Adipocytes',] %>% ggplot(percentage_concat_comnames, mapping=aes(x=Tissue, y=value, fill=Tissue))+ #fill=sex or fBMI
  geom_boxplot() + geom_boxplot(aes(color=Tissue),fatten = NULL, fill = NA) + theme_classic()+ scale_y_continuous(limits=c(0,70)) +
  scale_color_manual(values=c("#d95f0e","#fec44f"))+ scale_fill_manual(values=c("#d95f0e","#fec44f"))+
  ylab("Percentage %")+ theme_classic() +
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=18)) 

jpeg(paste0(output_path, "adipocytes_boxplot.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
print(p1)
dev.off()

# only 5 paired samples:
p1 <- paired_adipocytes[paired_adipocytes$cluster=='Adipocytes',] %>% ggplot(paired_adipocytes, mapping=aes(x=Tissue, y=value, fill=Tissue))+ #fill=sex or fBMI
  geom_boxplot() + geom_boxplot(aes(color=Tissue),fatten = NULL, fill = NA) + theme_classic()+ scale_y_continuous(limits=c(0,70)) +
  scale_color_manual(values=c("#d95f0e","#fec44f"))+ scale_fill_manual(values=c("#d95f0e","#fec44f"))+
  ylab("Percentage %")+
  theme_classic() +
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=18)) 

jpeg(paste0(output_path, "paired_adipocytes_boxplot.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# adipocytes: male and female
p1 <- percentage_concat_comnames[percentage_concat_comnames$cluster=='Adipocytes',] %>% ggplot(percentage_concat_comnames, mapping=aes(x=Tissue, y=value, fill=Sex))+ #fill=sex or fBMI
  geom_boxplot(width=0.75) + geom_boxplot(aes(color=Sex),fatten = NULL, fill = NA, width=0.75) + theme_classic()+ scale_y_continuous(limits=c(0,70)) +
  scale_color_manual(values=sex_colors)+ scale_fill_manual(values=sex_colors)+
  ylab("Percentage %") + theme_classic() +
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=18)) 

jpeg(paste0(output_path, "adipocytes_boxplot_F-M.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# only 5 paired samples:
p1 <- paired_adipocytes[paired_adipocytes$cluster=='Adipocytes',] %>% ggplot(paired_adipocytes, mapping=aes(x=Tissue, y=value, fill=Sex))+ #fill=sex or fBMI
  geom_boxplot() + geom_boxplot(aes(color=Sex),fatten = NULL, fill = NA) + theme_classic()+ scale_y_continuous(limits=c(0,70)) +
  scale_color_manual(values=sex_colors)+ scale_fill_manual(values=sex_colors)+
  ylab("Percentage %") + theme_classic() +
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=18)) 

jpeg(paste0(output_path, "paired_adipocytes_boxplot_F-M.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# adipocytes; BMI group
p1 <- percentage_concat_comnames[percentage_concat_comnames$cluster=='Adipocytes',] %>% ggplot(percentage_concat_comnames, mapping=aes(x=Tissue, y=value, fill=fBMI))+ #fill=sex or fBMI
  geom_boxplot(position = position_dodge(preserve = 'single')) + geom_boxplot(aes(color=fBMI),fatten = NULL, fill = NA, position = position_dodge(preserve = 'single')) + 
  theme_classic()+ scale_y_continuous(limits=c(0,70)) +
  scale_color_manual(values=BMI_colors)+ scale_fill_manual(values=BMI_colors)+
  ylab("Percentage %") + theme_classic() +
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=18)) 

jpeg(paste0(output_path, "adipocytes_boxplot_fBMI.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()
# only 5 paired samples:
p1 <- paired_adipocytes[paired_adipocytes$cluster=='Adipocytes',] %>% ggplot(paired_adipocytes, mapping=aes(x=Tissue, y=value, fill=fBMI))+ #fill=sex or fBMI
  geom_boxplot(position = position_dodge(preserve = 'single')) + geom_boxplot(aes(color=fBMI),fatten = NULL, fill = NA,position = position_dodge(preserve = 'single')) + 
  theme_classic()+ scale_y_continuous(limits=c(0,70)) +
  scale_color_manual(values=BMI_colors)+ scale_fill_manual(values=BMI_colors)+
  ylab("Percentage %") + theme_classic() +
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=18)) 

jpeg(paste0(output_path, "paired_adipocytes_boxplot_fBMI.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()

# adipocytes: age groups
p1 <- percentage_concat_comnames[percentage_concat_comnames$cluster=='Adipocytes',] %>% ggplot(percentage_concat_comnames, mapping=aes(x=Tissue, y=value, fill=fage))+ #fill=sex or fBMI
  geom_boxplot(position = position_dodge(preserve = 'single')) + geom_boxplot(aes(color=fage),fatten = NULL, fill = NA,position = position_dodge(preserve = 'single')) + 
  theme_classic()+ scale_y_continuous(limits=c(0,70)) +
  scale_color_manual(values=age_colors)+ scale_fill_manual(values=age_colors)+
  ylab("Percentage %") + theme_classic() +
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=18)) 

jpeg(paste0(output_path, "adipocytes_boxplot_fage.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()
# only 5 paired samples:
p1 <- paired_adipocytes[paired_adipocytes$cluster=='Adipocytes',] %>% ggplot(paired_adipocytes, mapping=aes(x=Tissue, y=value, fill=fage))+ #fill=sex or fBMI
  geom_boxplot(position = position_dodge(preserve = 'single')) + geom_boxplot(aes(color=fage),fatten = NULL, fill = NA,position = position_dodge(preserve = 'single')) + 
  theme_classic()+ scale_y_continuous(limits=c(0,70)) +
  scale_color_manual(values=age_colors)+ scale_fill_manual(values=age_colors)+
  ylab("Percentage %") + theme_classic() +
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size=18)) 

jpeg(paste0(output_path, "paired_adipocytes_boxplot_fage.jpeg"), width = 1500, height = 1500, quality = 100, bg = "transparent",res=300)
print(p1)
dev.off()




#### FIGURE 2:----
output_path <- "C:/Users/or711/OneDrive - post.bgu.ac.il/MSc-PhD/adipo/figs4paper/final_figs/fig2/"

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

