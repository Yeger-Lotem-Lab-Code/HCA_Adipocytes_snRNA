# Human subcutaneous and visceral adipocyte atlases uncover classical and non-classical adipocytes and depot-specific patterns

This repository contains data, code, and analysis as described in the paper ״Human subcutaneous and visceral adipocyte atlases uncover classical and non-classical adipocytes and depot-specific patterns".

## Table of contents
* [Article results](#article-results) 
    * [Pipeline code](#pipeline-code)
    * [Plots](#plots)
* [Data availability](#Data-availability)

### Article results
Results are found in their respective "article_results/Figure X" folder.

#### Pipeline code
the main code was written in python. It includes data processing, machine learning construction and validation. 
The pipeline code can be reviewed and used by running main.py from src folder.

#### Plots
Data to generate each panel is presented in "Data source file.xlsx" and in their respective "article_results/Figure X" folder.
Plots were generated as part of the pipeline via R.

| Figure    | Panel | Scripts                                                                                                                   |
|-----------|-------|----------------------------------------------------------------------------------------------------------------------------|
| Main_figs_info |       | Main_figs.info.R (includes necessary libraries, Seurat objects)                                                       |
| Figure1 - Single-nuclei atlases of human subcutaneous (hSAT) and visceral (hVAT) adipose tissues | A-D   | Fig1.R                      |
| Figure2   | All   | Fig2.R                                                                                                                     |
| Figure3   | All   | GSEA_enrichment.R → gsea_color_new.py → Fig3.R                                                                             |
| Figure4   | All   | Fig4.R                                                                                                                                                                                                                          |
| Figure5   | All     | Fig5.R                                                                                         |
|    Figure6       | A-C, E-G  | Fig6.R                                                                                                             |
|   | D,H  | palantir_analysis.py                                                                                                    |
| Figure7   | A     | cellphone_pipeline.py                                                   |
|           | B-I   | cellchat.R                                                                                       |

    


#### Data availability
The dataset is available on Zenode - (https://cellxgene.cziscience.com/collections/ba84c7ba-8d8c-4720-a76e-3ee37dc89f0b)


Analysis pipeline files according to order of usage:
1. raw_objects.R - create Seurat objects per sample after CellBender AND CellRanger filtration for ambient RNA, and nFeature (=>200) and mitoRNA cutoff (=<20%). 
2. basic_analysis.R - perform clustering, annotation of adipocyte clusters and create violin plots for nCounts and nFeatures per cluster.
3. remove_lowQuality.R - removes clusters with mean nFeatures (nCounts) in cluster <= (mean nFeatures (nCounts) - S.D.) in a sample.
4. doubletFinder.R - doubletFinder pipeline analysis for each Seurat object.
5. harmony_integration.R - integration analysis with harmony.
6. AssessNodes.R - Seurat V2 function repurposed to fit Seurat V3 objects.
7. GSEA_code.R - run Gene Set Enrichment Analysis with clusterProfiler for each subset in Seurat object.
