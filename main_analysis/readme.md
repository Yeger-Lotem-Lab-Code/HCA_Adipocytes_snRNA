

# Analysis pipeline files according to order of usage:
1. raw_objects.R - create Seurat objects per sample after CellBender AND CellRanger filtration for ambient RNA, and nFeature (=>200) and mitoRNA cutoff (=<20%). 
2. basic_analysis.R - perform clustering, annotation of adipocyte clusters and create violin plots for nCounts and nFeatures per cluster.
3. remove_lowQuality.R - removes clusters with mean nFeatures (nCounts) in cluster <= (mean nFeatures (nCounts) - S.D.) in a sample.
4. doubletFinder.R - doubletFinder pipeline analysis for each Seurat object.
5. harmony_integration.R - integration analysis with harmony.
6. AssessNodes.R - Seurat V2 function repurposed to fit Seurat V3 objects.
7. GSEA_code.R - run Gene Set Enrichment Analysis with clusterProfiler for each subset in Seurat object.
