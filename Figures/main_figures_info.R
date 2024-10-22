# scripts for each main fig:

# libraries:
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(dplyr))
library(plyr)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggridges))
library(patchwork)
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(hdf5r))
# suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(UniprotR))
library(rstatix)
library(stats)

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
cols4bars <- c("#00441b","#006d2c","#238b45","#41ab5d","#74c476","#a1d99b","#c7e9c0",
               "#49006a","#7a0177","#ae017e","#dd3497","#f768a1","#fa9fb5","#fcc5c0", "#fde0dd")
cols4bars_vis_emont2 <- c("#969696",
                          "#49006a","#7a0177","#ae017e","#dd3497","#f768a1","#fa9fb5","#fcc5c0", "#fde0dd")
cols4bars_sc_emont2 <- c("#969696",
                         "#00441b","#006d2c","#238b45","#41ab5d","#74c476","#a1d99b","#c7e9c0")
sc_aspc_adipo_col <- c("SA1" = "#ff7f00", "SA2" = "#6a3d9a", "SA3" = "#6a3d9a","SA4" = "#6a3d9a", 
                   "SA5" = "#6a3d9a", "SA6" = "#6a3d9a", "SA7" = "#6a3d9a",
                   "ASPC 1"= "#33a02c", "ASPC 2" = "#33a02c", "ASPC 3"= "#33a02c", "ASPC 4" = "#33a02c",  "ASPC 5" = "#33a02c")

vis_aspc_adipo_col <- c("VA1" = "#ff7f00", "VA2" = "#6a3d9a", "VA3" = "#6a3d9a",
                        "VA4" = "#6a3d9a", "VA5" = "#6a3d9a", "VA6" = "#6a3d9a", 
                        "VA7" = "#6a3d9a","VA8" = "#6a3d9a",
                        "ASPC 1"= "#33a02c", "ASPC 2" = "#33a02c", "ASPC 3"= "#33a02c", "ASPC 4" = "#33a02c")

aspc_only <- c("ASPC 1"= "#1f78b4", "ASPC 2" = "#33a02c", "ASPC 3"= "#a6cee3", "ASPC 4" = "#b2df8a",  "ASPC 5" = "#fb9a99")

#### seurat objects:----
obj_path <- "/seurat/objects/"
vis10 <- readRDS(paste0(obj_path, "vis10_diet.rds"))
sc5 <- readRDS(paste0(obj_path, "sc5_diet.rds"))
adipo10 <- readRDS(paste0(obj_path, "adipo10_diet.rds"))
adipo5 <- readRDS(paste0(obj_path, "adipo5_diet.rds"))
adipo_vissc <- readRDS(paste0(obj_path, "adipo_vissc_diet.rds")) #both VA1-8 and SA1-7 are under adipo_ann
vis_house_emont <- readRDS(paste0(obj_path, "vis_house_emont.rds"))
sc_house_emont <- readRDS(paste0(obj_path, "sc_house_emont.rds"))
aspc_adipo10 <- readRDS(paste0(obj_path, "vis_ASPC_adipo.rds"))
Idents(object = aspc_adipo10) <- "final_detailed"
aspc_vis_only <- subset(aspc_adipo10, idents = c("ASPC 1", "ASPC 2", "ASPC 3", "ASPC 4"), invert=F)

aspc_adipo5 <- readRDS(paste0(obj_path, "sc_ASPC_adipo.rds"))
Idents(object = aspc_adipo5) <- "final_detailed"
aspc_sc_only <- subset(aspc_adipo5, idents = c("ASPC 1", "ASPC 2", "ASPC 3", "ASPC 4", "ASPC 5"), invert=F)

#### table preparations:----
# in-house meta-data:
tables_path <- "tables/path/"
house_meta <-  read.csv(paste0(tables_path, "clinical_parameters_2024-10-20.csv"), header=TRUE)
house_meta$Age <- round(house_meta$Age, digits = 1)

