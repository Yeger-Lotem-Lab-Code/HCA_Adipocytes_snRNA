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
Data to generate each panel is deposited on Cellxgene if it requires the seurat/cellchat objects, and relevant CSV files are inside the figure's folder.
Plots were generated as part of the pipeline via R.

| Figure    | Panel    | Scripts                                                                                  |
|-----------|----------|-------------------------------------------------------------------------------------------|
| Main_figs_info |          | Main_figs.info.R (includes necessary libraries, Seurat objects)                      |
| Figure1 - Single-nuclei atlases of human subcutaneous (hSAT) and visceral (hVAT) adipose tissues | A-D | Fig1.R |
| Figure2   | All      | Fig2.R                                                                                   |
| Figure3   | All      | GSEA_enrichment.R → gsea_color_new.py → Fig3.R                                            |
| Figure4   | All      | Fig4.R                                                                                    |
| Figure5   | All      | Fig5.R                                                                                    |
| Figure6   | A-C, E-G | Fig6.R                                                                                    |
|           | D, H     | palantir_analysis.py                                                                      |
| Figure7   | A        | cellphone_pipeline.py                                                                     |
|           | B-I      | cellchat.R                                                                                |
                                                                       |

    


#### Data availability
The dataset is available on Cellxgene - (https://cellxgene.cziscience.com/collections/ba84c7ba-8d8c-4720-a76e-3ee37dc89f0b)

