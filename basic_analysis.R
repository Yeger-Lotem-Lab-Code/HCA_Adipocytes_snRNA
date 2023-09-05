library(Seurat)
library(dplyr)
library(ggplot2)
library(plotly)
library(ggpubr)

output_path = ".../basic_analysis/"
load_path = ".../raw_objects/"

basic_analysis <- function(sample_name){
  dir.create(paste(output_path, sample_name, sep = ""))
  obj_seurat <- readRDS(paste(load_path, sample_name,"/", sample_name,".rds", sep=""))
  obj_seurat <- NormalizeData(obj_seurat)
  obj_seurat <- FindVariableFeatures(obj_seurat, nfeatures = 2000)
  obj_seurat <- ScaleData(obj_seurat, features = rownames(obj_seurat))
  obj_seurat <- RunPCA(obj_seurat)
  obj_seurat <- JackStraw(obj_seurat, dims = 50)
  obj_seurat <- ScoreJackStraw(obj_seurat, dims = 1:50, reduction = "pca")
  plt1 <- JackStrawPlot(obj_seurat, dims = 1:50)
  jpeg(paste(output_path, sample_name,"/", sample_name, "_jackstrawPlot.jpeg", sep = ""), 
       width = 1000, height = 750)
  print(plt1)
  dev.off()
  rm(plt1)
  saveRDS(obj_seurat, file = (paste(output_path, sample_name,"/", sample_name,".rds", sep="")))
  rm(obj_seurat)
}

cluster_UMAP <- function(sample_name, num_dims){
  obj_seurat <- readRDS(paste(output_path,"/", sample_name,"/", sample_name,".rds", sep=""))
  obj_seurat <- FindNeighbors(obj_seurat, dims = 1:num_dims)
  obj_seurat <- FindClusters(obj_seurat, resolution = 0.8)
  obj_seurat <- RunUMAP(obj_seurat, dims = 1:num_dims)
  
  p1 <- FeaturePlot(obj_seurat, features = "nCount_RNA", pt.size = 1.5)
  p2 <- FeaturePlot(obj_seurat, features = "nFeature_RNA", pt.size = 1.5)
  p3 <- FeaturePlot(obj_seurat, features = "percent.mt", pt.size = 1.5)
  p4 <- FeaturePlot(obj_seurat, features = "AQP7", pt.size = 1.5)
  jpeg(paste(output_path, sample_name,"/", sample_name, "_4featurePlots.jpeg", sep = ""), 
       width = 3000, height = 3000)
  print(ggarrange(p1,p2,p3,p4, ncol=2, nrow=2))
  dev.off()
  rm(p1,p2,p3,p4)
  
  plt2 <- DimPlot(obj_seurat, reduction = "umap", label = T, pt.size = 1.5, label.size = 8)+
    theme(legend.position = "none")
  jpeg(paste(output_path, sample_name,"/", sample_name, "_umapByCluster.jpeg", sep = ""), 
       width = 1500, height = 1250)
  print(plt2)
  dev.off()
  rm(plt2)
  
  plt3 <- VlnPlot(obj_seurat, features = "nCount_RNA", pt.size = 0)+ 
    geom_boxplot(width=0.3)+
    theme(legend.position = 'none') +
    scale_y_continuous(trans = "log10")
  plt4 <- VlnPlot(obj_seurat, features = "nFeature_RNA", pt.size = 0)+ 
    geom_boxplot(width=0.3)+
    theme(legend.position = 'none') +
    scale_y_continuous(trans = "log10")
  plt5 <- VlnPlot(obj_seurat, features = "percent.mt", pt.size = 0)+ 
    geom_boxplot(width=0.3)+
    theme(legend.position = 'none') +
  jpeg(paste(output_path, sample_name,"/", sample_name, "_violinPlots.jpeg", sep = ""), 
       width = 1000, height = 3000)
  print(plt3+plt4+plt5)
  dev.off()
  rm(plt4, plt3,plt5)
  
  saveRDS(obj_seurat, file = (paste(output_path, sample_name,"/", sample_name,".rds", sep="")))
  rm(obj_seurat)
}

#
sample_list <- list.dirs(load_path, full.names = F, recursive = F)
for (i in 3:length(sample_list)) {
  print(sample_list[i])
  basic_analysis(sample_name = sample_list[i])
}

#
dim_list <- c() # enter dim number for each sample according to jackstrawplot
for (i in 3:length(sample_list)) {
  print(sample_list[i])
  cluster_UMAP(sample_name = sample_list[i], num_dims = dim_list[i-2])
}
