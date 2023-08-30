suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(DoubletFinder))

output_path = "/gpfs0/estiyl/users/mziv/immune/samples_29.1/remove_low_quality/"
input_path = "/gpfs0/estiyl/users/mziv/immune/samples_29.1/basic_analysis/"

summary_DF <- data.frame(row.names = c("num cells before",
                                       "max nFeature", "mean nFeature", "min nFeature",
                                       "max nCounts","mean nCounts", "min nCounts",
                                       "max percent.mt","mean percent.mt", "min percent.mt", "clusters unfit",
                                       "num cells after"))


sample_list <- list.dirs(input_path, full.names = F, recursive = F)
for (i in 1:length(sample_list)) {
  print(sample_list[i])
  dir.create(paste(output_path, "/", sample_list[i], sep = ""))
  seurat_obj <- readRDS(file = paste(input_path, sample_list[i],"_cluster.rds", sep = ""))
  col2add <- list(dim(seurat_obj)[2], #number of cells
                  max(seurat_obj@meta.data$nFeature_RNA), mean(seurat_obj@meta.data$nFeature_RNA), min(seurat_obj@meta.data$nFeature_RNA),
                  max(seurat_obj@meta.data$nCount_RNA), mean(seurat_obj@meta.data$nCount_RNA), min(seurat_obj@meta.data$nCount_RNA),
                  max(seurat_obj@meta.data$percent.mt), mean(seurat_obj@meta.data$percent.mt), min(seurat_obj@meta.data$percent.mt))
  temp_df <- data.frame(row.names = levels(seurat_obj@active.ident))
  #for loop to find the mean features and counts per cluster
  for(j in row.names(temp_df)) {
    temp_df[j,"mean_counts"] <-  mean(seurat_obj@meta.data$nCount_RNA[seurat_obj@meta.data$seurat_clusters==j])
    temp_df[j,"mean_Features"] <-  mean(seurat_obj@meta.data$nFeature_RNA[seurat_obj@meta.data$seurat_clusters==j])
  }
  # calculate threshold for nFeatures and nCounts
  features_thresh <- (mean(temp_df$mean_Features) - sd(temp_df$mean_Features))
  counts_thresh <- (mean(temp_df$mean_counts) - sd(temp_df$mean_counts))
  # decide which clusters are unfit and remove them
  toRemove <- list()
  for(k in row.names(temp_df)){
    if(temp_df[k, "mean_Features"]< features_thresh & temp_df[k, "mean_counts"]< counts_thresh){
      toRemove <- c(toRemove, k) # list of clusters with nFeatures or nCounts lower than the threshold
    }
  }
  # #remove clusters with low expression
  seurat_obj <- subset(seurat_obj, idents = toRemove, invert = T)
  toRemove <- paste(toRemove, collapse=', ' ) # colapse to create one 'string' of cluster numbers
  # check for #cells after removal
  col2add <- c(col2add, toRemove,dim(seurat_obj)[2])
  # insert relevant data to summary table
  summary_DF[[sample_list[i]]] <- col2add
  saveRDS(seurat_obj, file = paste(output_path, sample_list[i],"/", sample_list[i], ".rds", sep = ""))
  # rm(j,k,toRemove, temp_df, col2add, seurat_obj)
  rm(j,k,toRemove, temp_df, col2add, seurat_obj)

}
summary_DF <- apply(summary_DF, 2, as.character)
row.names(summary_DF) <- c("num cells before",
                           "max nFeature", "mean nFeature", "min nFeature",
                           "max nCounts","mean nCounts", "min nCounts",
                           "max percent.mt","mean percent.mt", "min percent.mt","clusters unfit",
                           "num cells after")
write.csv(summary_DF, file = paste(output_path, "summary_DF.csv", sep = ""), row.names = TRUE)



# basic_analysis <- function(sample_name){
#   #dir.create(paste(output_path, sample_name, sep = ""))
#   obj_seurat <- readRDS(paste(output_path, sample_name,".rds", sep=""))
#   obj_seurat[["RNA_snn_res.0.8"]] <- NULL ####
#   obj_seurat[["seurat_clusters"]] <- NULL ####
#   obj_seurat <- NormalizeData(obj_seurat)
#   obj_seurat <- FindVariableFeatures(obj_seurat, nfeatures = 2000)
#   obj_seurat <- ScaleData(obj_seurat, features = rownames(obj_seurat))
#   obj_seurat <- RunPCA(obj_seurat)
#   obj_seurat <- JackStraw(obj_seurat, dims = 50)
#   obj_seurat <- ScoreJackStraw(obj_seurat, dims = 1:50, reduction = "pca")
#   plt1 <- JackStrawPlot(obj_seurat, dims = 1:50)
#   jpeg(paste(output_path, sample_name, "_jackstrawPlot.jpeg", sep = ""),
#        width = 1000, height = 750)
#   print(plt1)
#   dev.off()
#   rm(plt1)
#   saveRDS(obj_seurat, file = (paste(output_path, sample_name,".rds", sep="")))
#   rm(obj_seurat)
# }
# 
# cluster_UMAP <- function(sample_name, num_dims){
#   obj_seurat <- readRDS(paste(output_path, sample_name,".rds", sep=""))
#   obj_seurat <- FindNeighbors(obj_seurat, dims = 1:num_dims)
#   obj_seurat <- FindClusters(obj_seurat, resolution = 0.8)
#   obj_seurat <- RunUMAP(obj_seurat, dims = 1:num_dims)
# 
#   p1 <- FeaturePlot(obj_seurat, features = "nCount_RNA", pt.size = 1.5)
#   p2 <- FeaturePlot(obj_seurat, features = "nFeature_RNA", pt.size = 1.5)
#   p3 <- FeaturePlot(obj_seurat, features = "percent.mt", pt.size = 1.5)
#   jpeg(paste(output_path, sample_name, "_3featurePlots.jpeg", sep = ""),
#        width = 3000, height = 3000)
#   print(ggarrange(p1,p2,p3, ncol=2, nrow=2))
#   dev.off()
#   rm(p1,p2,p3)
# 
#   plt2 <- DimPlot(obj_seurat, reduction = "umap", label = T, pt.size = 1.5, label.size = 8)+
#     theme(legend.position = "none")
#   jpeg(paste(output_path, sample_name, "_umapByCluster.jpeg", sep = ""),
#        width = 1500, height = 1250)
#   print(plt2)
#   dev.off()
#   rm(plt2)
# 
#   plt3 <- VlnPlot(obj_seurat, features = "nCount_RNA", pt.size = 0)+
#     geom_boxplot(width=0.3)+
#     theme(legend.position = 'none') +
#     scale_y_continuous(trans = "log10")
#   plt4 <- VlnPlot(obj_seurat, features = "nFeature_RNA", pt.size = 0)+
#     geom_boxplot(width=0.3)+
#     theme(legend.position = 'none') +
#     scale_y_continuous(trans = "log10")
#   plt5 <- VlnPlot(obj_seurat, features = "percent.mt", pt.size = 0)+
#     geom_boxplot(width=0.3)+
#     theme(legend.position = 'none') +
#     scale_y_continuous(trans = "log10")
# 
#   jpeg(paste(output_path, sample_name, "_4violinPlots.jpeg", sep = ""),
#        width = 3000, height = 3000)
#   print(ggarrange(plt3,plt4,plt5, ncol=2, nrow=2))
#   dev.off()
#   rm(plt4, plt3,plt5)
# 
#   saveRDS(obj_seurat, file = (paste(output_path, sample_name,".rds", sep="")))
#   rm(obj_seurat)
# }

# 
# 
# basic_analysis(sample_name = sample_name)



#cluster_UMAP(sample_name = sample_name, num_dims = 34)

