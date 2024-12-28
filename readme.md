# Human subcutaneous and visceral adipocyte atlases uncover classical and non-classical adipocytes and depot-specific patterns

This repository contains code for the analyses described in the paper ״Human subcutaneous and visceral adipocyte atlases uncover classical and non-classical adipocytes and depot-specific patterns".

## Table of contents
* [Article results](#article-results) 
    * [Pipeline code](#pipeline-code)
    * [Figures](#Figures)
* [Data availability](#Data-availability)

### Article results
Results are found in their respective "Figure X" folder.

#### Pipeline code
The main code was written in R. It includes data processing, quality control, clustering and annotating the code file are under the "main_analysis" folder with the description of the files. 

#### Figures
The data used to generate each panel is available on CELLxGENE, and the corresponding source files can be found within the figure's folder.

# Figures and Scripts Overview

| **Figure**               | **Panel**    | **Scripts**                                                                                  |
|---------------------------|--------------|---------------------------------------------------------------------------------------------|
| **Main_figs_info**        |              | `Main_figs.info.R` (includes necessary libraries, Seurat objects)                           |
| **Figure1**               | A-D          | `Fig1.R`                                                                                     |
| **Figure2**               | All          | `Fig2.R`                                                                                   |
| **Figure3**               | All          | `GSEA_enrichment.R` → `gsea_color_new.py` → `Fig3.R`                                       |
| **Figure4**               | All          | `Fig4.R`                                                                                   |
| **Figure5**               | All          | `Fig5.R`                                                                                   |
| **Figure6**               | A-C, E-G     | `Fig6.R`                                                                                   |
|                           | D, H         | `palantir_analysis.py`                                                                      |
| **Figure7**               | A            | `cellphone_pipeline.py`                                                                     |
|                           | B-I          | `cellchat.R`                                                                                |
| **Ex_Data_Fig.1**         | All          | `Ex_Data_fig1.R`                                                                            |
| **Ex_Data_Fig.4**         | All          | `Ex_Data_fig4.R`                                                                            |
| **Ex_Data_Fig.6**         | All          | `Ex_Data_fig6.R`                                                                            |
| **Ex_Data_Fig.7**         | All          | `Ex_Data_fig7.R`                                                                            |
| **Supp_Fig1**             | All          | `Fig.S1.R`                                                                                 |
| **Supp_Fig2**             | All          | `Fig.S2.R`                                                                                 |
| **Supp_Fig3+4**           | All          | `Fig.S3+S4.R`                                                                              |
| **Supp_Fig5+6+7**         | All          | `Fig.S5+S6+S7.R`                                                                           |
| **Supp_FigS8**            | All          | `river_contigency.R`                                                                        |
| **Supp_FigS9**            | All          | `Fig.S9.R`                                                                                 |
| **Supp_FigS10**           | All          | `spatial.R`                                                                                |
| **Supp_FigS11**           | All          | `Fig.S11.R`                                                                                |
| **Supp_FigS12-14**        | All          | `Fig.S12-S14.R`                                                                            |
| **Supp_FigS15**           | All          | `Fig.S15.R`                                                                                |
| **Supp_FigS116-17**       | All          | `Fig.S116-17.R`                                                                            |
    


#### Data availability
The dataset is available on CELLxGENE - (https://cellxgene.cziscience.com/collections/ba84c7ba-8d8c-4720-a76e-3ee37dc89f0b)

## Our repository is also on Zenodo : (https://zenodo.org/records/14001210)
