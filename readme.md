# Human subcutaneous and visceral adipocyte atlases uncover classical and non-classical adipocytes and depot-specific patterns

This repository contains data, code, and analysis as described in the paper ״Human subcutaneous and visceral adipocyte atlases uncover classical and non-classical adipocytes and depot-specific patterns".
Link to data: (https://cellxgene.cziscience.com/collections/ba84c7ba-8d8c-4720-a76e-3ee37dc89f0b)
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

| Figure    | Panel | Script                                                                                                                   |
|---------- |-------|------------------------------------------------------------------------------------------------------------------------- |
| Figure1   | A     | build_aneuploidy_table_GISTIC.py (pipeline)                                                                              |
| Figure2   | A-D   | shap_analysis.py (pipline)                                                                                               |
| Figure2   | E     | correlation_analysis.py (pipline)                                                                                        |
| Figure3   | A-E   | CCL_main.py (pipline)                                                                                                    |
| Figure3   | F     | correlation_analysis.py (pipline)                                                                                        |
| Figure4   | B     | Fig4B_13q_gain_prevalence_barplot                                                                                        |
|           | D,F   | violin_plots.R                                                                                                           |
|           | H,J   | Fig4H-J_boxplots                                                                                                         |
| Figure5   | A     | correlation_analysis.py (pipline)                                                                                        |
| Figure5   | D,E   | violin_plots.R                                                                                                           |
|           | G,H   | Fig5G-H_cooccurrence.R                                                                                                   |
| FigureS3  | A-C   | Correlation analysis: article_results\Figure S3\correlation analysis.py                                                  |
|           |       | Plots: article_results\Figure S3\Fig_S3_plots.r                                                                          |
| FigureS4  | A-B   | curves - compare_models.py (pipline)                                                                                     |
|           |       | boxplot - article_results\Figure S4\models_performance.R                                                                 |
| FigureS5  | A,B   | shap_analysis.py (pipeline)                                                                                              |
| FigureS6  |       | correlation_analysis.py (pipline)                                                                                        |   
| FigureS7  |       | construct dataset with other threshold: src\other_threshold_percentage.py                                                |  
|           |       | To find best model and to plot SHAP, run compare_models.py and SHAP analysis respectively with the corresponding dataset | 
| FigureS8  | A     | parameter_tunning.py (pipline)                                                                                           |
| FigureS8  | B-C   | For model performance and SHAP analysis, run compare_models.py and SHAP analysis with best hyperpatameter per model      |
| FigureS9  | A     | curves - compare_models.py (pipline)                                                                                     |
|           |       | boxplot - article_results\Figure S4\models_performance.R                                                                 |
| FigureS9  | B-C   | shap_analysis.py (pipeline)                                                                                              |
| FigureS10 | A-B   | article_results\Figure S10\top_driver_analysis.py                                                                        |
| FigureS11 | A     | To find best model and to plot SHAP, run compare_models.py and SHAP analysis respectively with the corresponding dataset | 
| FigureS11 | B-C   | article_results\Figure S11\whole_chromosome_main.py                                                                      |
| FigureS12 |       | article_results\src\boxplots.R                                                                                           |
| FigureS13 | A     | article_results\src\bar_plots.R                                                                                          |
| FigureS14 | C     | article_results\src\boxplots.R                                                                                           |
| FigureS15 | A     | article_results\src\bar_plots.R                                                                                          |
|           | D     | article_results\src\boxplots.R                                                                                           |
| FigureS16 | A     | article_results\src\bar_plots.R                                                                                          |
|           | D     | article_results\src\boxplots.R                                                                                           |
| FigureS17 | A-B   | paralogs_analysis.py (pipline)         


#### Data availability
The dataset is available on Zenode - https://zenodo.org/record/8199048
Analysis pipeline files according to order of usage:
1. raw_objects.R - create Seurat objects per sample after CellBender AND CellRanger filtration for ambient RNA, and nFeature (=>200) and mitoRNA cutoff (=<20%). 
2. basic_analysis.R - perform clustering, annotation of adipocyte clusters and create violin plots for nCounts and nFeatures per cluster.
3. remove_lowQuality.R - removes clusters with mean nFeatures (nCounts) in cluster <= (mean nFeatures (nCounts) - S.D.) in a sample.
4. doubletFinder.R - doubletFinder pipeline analysis for each Seurat object.
5. harmony_integration.R - integration analysis with harmony.
6. AssessNodes.R - Seurat V2 function repurposed to fit Seurat V3 objects.
7. GSEA_code.R - run Gene Set Enrichment Analysis with clusterProfiler for each subset in Seurat object.
