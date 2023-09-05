library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(plotly)
library(ggpubr)

input_path = ".../lowQuality_removed/"
output_path = "/doubletFinder_analysis/"

### ------------ functions ------------- ###
basic_analysis <- function(sample_name, list2remove){
  obj_seurat <- readRDS(paste(input_path, sample_name,"/", sample_name,".rds", sep=""))
  obj_seurat <- subset(obj_seurat, idents = list2remove, invert = T)
  obj_seurat[["RNA_snn_res.0.5"]] <- NULL
  obj_seurat[["seurat_clusters"]] <- NULL
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
  obj_seurat <- FindClusters(obj_seurat, resolution = 0.5)
  obj_seurat <- RunUMAP(obj_seurat, dims = 1:num_dims)
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
  plt5 <- VlnPlot(obj_seurat, features = c("AQP7", "ADIPOQ"))
  jpeg(paste(output_path, sample_name,"/", sample_name, "_violinPlots.jpeg", sep = ""), 
       width = 1000, height = 1800)
  print(plt3+plt4)
  dev.off()
  jpeg(paste(output_path, sample_name,"/", sample_name, "_vlnPlot_adipocytes.jpeg", sep = ""), 
       width = 2000, height = 1000)
  print(plt5)
  dev.off()
  rm(plt5, plt4, plt3)
  saveRDS(obj_seurat, file = (paste(output_path,sample_name,"/", sample_name,".rds", sep="")))
  rm(obj_seurat)
}

find_pk <- function(sample_name, num_dims){
  obj_seurat <- readRDS(paste(input_path, sample_name,"/", sample_name,".rds", sep=""))
  sweep.res.list <- paramSweep_v3(obj_seurat, PCs = 1:num_dims, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric),]$pK))
  nExp_poi <- round(0.075*nrow(obj_seurat@meta.data))
  annotations <- obj_seurat@active.ident ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) ## Assuming 7.5% doublet formation rate:
  list2return <- list("pK" = pK, "nExp_poi" = nExp_poi, "nExp_poi.adj" = nExp_poi.adj)
  rm(sweep.res.list, sweep.stats, bcmvn, pK, nExp_poi, annotations, homotypic.prop, nExp_poi.adj, obj_seurat)
  return(list2return)
}

# lists-------
sample_list <- list.dirs(input_path, full.names = F, recursive = F)

# dims after removing low-quality:
dim_list <- list()

DB_params <- data.frame(row.names = c("pK", "nExp_poi", "nExp_poi.adj"))

# loops-------
# run after removing low-quality clusters
for (i in 1:length(sample_list)) {
  print(sample_list[i])
  DB_params[[sample_list[[i]]]] <- find_pk(sample_name = sample_list[i], num_dims = dim_list[[i]])
}


DB_params <- data.frame(lapply(DB_params, as.character), stringsAsFactors=FALSE)

for (i in 1:length(sample_list)) {
  print(sample_list[i])
  dir.create(paste0(output_path, sample_list[i], "/"))
  obj_seurat <- readRDS(paste(input_path, sample_list[i],"/", sample_list[i],".rds", sep=""))
  obj_seurat <- doubletFinder_v3(obj_seurat, PCs = 1:(dim_list[[i]]), pN = 0.25, pK = DB_params[[1,i]], 
                                 nExp = DB_params[[2,i]], reuse.pANN = FALSE, sct = FALSE)
  obj_seurat <- doubletFinder_v3(obj_seurat, PCs = 1:(dim_list[[i]]), pN = 0.25, pK = DB_params[[1,i]], nExp = DB_params[[3,i]], 
                                 reuse.pANN = paste0("pANN_0.25_", as.character(DB_params[[1,i]]),"_", as.character(DB_params[[2,i]])), sct = FALSE)
  saveRDS(obj_seurat, file = (paste(output_path, sample_list[i],"/", sample_list[i],".rds", sep="")))
  rm(obj_seurat)
}

# list the cells with high doublet score
#1. load each object 
#2. write csv table with; cell names, score and DF classification
#3. write csv for clusters percentage of doublets
#4. figures of umap colored by clssification1 and classification2 and by score

for (i in 1:length(sample_list)) {
  print(sample_list[i])
  obj_seurat <- readRDS(paste(output_path, sample_list[i],"/", sample_list[i],".rds", sep=""))
  pANN <- paste0("pANN_0.25_", as.character(DB_params[[1,i]]),"_", as.character(DB_params[[2,i]]))
  DF.classification1 <- paste0("DF.classifications_0.25_", as.character(DB_params[[1,i]]),"_", as.character(DB_params[[2,i]]))
  DF.classification2 <- paste0("DF.classifications_0.25_", as.character(DB_params[[1,i]]),"_", as.character(DB_params[[3,i]]))
  DF_classification <- data.frame(row.names = row.names(obj_seurat@meta.data), obj_seurat@active.ident,
                                  obj_seurat@meta.data[,pANN],
                                  obj_seurat@meta.data[,DF.classification1],
                                  obj_seurat@meta.data[,DF.classification2])
  colnames(DF_classification) <- c("cluster", "pANN", "DF.classification", "DF.classification_homotypic")
  write.csv(DF_classification, file = paste(output_path, sample_list[i], "_DF_classifications.csv", sep = ""), row.names = TRUE)
  DF_summary <- table(DF_classification$cluster, DF_classification$DF.classification_homotypic)
  DF_summary <- prop.table(DF_summary, 1)*100
  write.csv(DF_summary, file = paste(output_path, sample_list[i], "/", sample_list[i], "_DF_clusters_percentage.csv", sep = ""), row.names = TRUE)
  plt1 <- FeaturePlot(obj_seurat, reduction = "umap", pt.size = 0.8, features = pANN)
  plt2 <- DimPlot(obj_seurat, reduction = "umap", pt.size = 0.8, group.by = DF.classification1, cols = c("purple", "grey"))
  plt3 <- DimPlot(obj_seurat, reduction = "umap", pt.size = 0.8, group.by = DF.classification2, cols = c("purple", "grey"))
  plt4 <- VlnPlot(obj_seurat, features = pANN)
  jpeg(paste(output_path, sample_list[i],"/", sample_list[i], "_DF_score_plots.jpeg", sep = ""),
       width = 3000, height = 3000)
  print(ggarrange(plt1,plt2,plt3,plt4, ncol=2, nrow=2))
  dev.off()
  rm(obj_seurat,DF_summary, pANN, DF_classification, DF.classification1, DF.classification2,plt1,plt2,plt3,plt4) 
}

# read both tables from last loop:
# remove all cells defined as doublets by DF_homotypic (DF_classification2)
# remove clusters with doublets>=65%
# clear object and save
# export csv table with summary

input_path = ".../doubletFinder_analysis/"
output_path = ".../doubletFinder_analysis/doublets_removed/"
summary_DF <- data.frame(row.names = c("num cells before", "clusters removed", "num cells removed", "num cells after"))

for (i in 1:length(sample_list)) {
  print(sample_list[i])
  seurat_obj <- readRDS(paste(input_path, sample_list[i],"/", sample_list[i],".rds", sep=""))
  DF_classifications <- read.csv(file = paste0(input_path, sample_list[i], "_DF_classifications.csv"), header = T, row.names = 1)
  clusters_percentage <- read.csv(file = paste0(input_path, sample_list[i], "/", sample_list[i], "_DF_clusters_percentage.csv"), header = T, row.names = 1)
  cells2remove <- row.names(DF_classifications[which(DF_classifications$DF.classification_homotypic=="Doublet"),])
  clusters2remove <- as.numeric(row.names(clusters_percentage[which(clusters_percentage$Doublet>=65),]))
  toRemove <- paste(clusters2remove, collapse=', ' )
  col2add <- list(dim(seurat_obj)[2], toRemove)
  seurat_obj <- subset(seurat_obj, cells = cells2remove, invert = T)
  if (length(clusters2remove)>0) {
    seurat_obj <- subset(seurat_obj, idents = clusters2remove, invert = T)
  }
  num_cells_removed <- c(col2add[[1]] - dim(seurat_obj)[2])
  col2add <- c(col2add, num_cells_removed, dim(seurat_obj)[2])
  summary_DF[[sample_list[[i]]]] <- col2add
  #clear seurat object and save:
  seurat_obj <- DietSeurat(seurat_obj, counts = T, data = T, scale.data = F)
  seurat_obj@meta.data$RNA_snn_res.0.8 <- NULL
  seurat_obj@meta.data$seurat_clusters <- NULL
  pANN <- paste0("pANN_0.25_", as.character(DB_params[[1,i]]),"_", as.character(DB_params[[2,i]]))
  DF.classification1 <- paste0("DF.classifications_0.25_", as.character(DB_params[[1,i]]),"_", as.character(DB_params[[2,i]]))
  DF.classification2 <- paste0("DF.classifications_0.25_", as.character(DB_params[[1,i]]),"_", as.character(DB_params[[3,i]]))
  seurat_obj@meta.data[,pANN] <- NULL
  seurat_obj@meta.data[,DF.classification1] <- NULL
  seurat_obj@meta.data[,DF.classification2] <- NULL
  saveRDS(seurat_obj, file = paste(output_path, sample_list[i], ".rds", sep = ""))
  rm(seurat_obj, DF_classifications, cells2remove, clusters2remove, toRemove, col2add, clusters_percentage,
     num_cells_removed,pANN,DF.classification1,DF.classification2)
}

summary_DF <- apply(summary_DF, 2, as.character)
row.names(summary_DF) <- c("num cells before", "clusters removed", "num cells removed", "num cells after")
write.csv(summary_DF, file = paste(output_path, "summary_DF.csv", sep = ""), row.names = TRUE)
