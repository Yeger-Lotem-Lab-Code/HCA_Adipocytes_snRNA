library(Seurat)
library(dplyr)
library(ggplot2)
library(plotly)

input_path = ".../basic_analysis/"
output_path = ".../lowQuality_removed/"

summary_DF <- data.frame(row.names = c("num cells before",
                                       "max nFeature", "mean nFeature", "min nFeature", 
                                       "max nCounts","mean nCounts", "min nCounts", 
                                       "max percent.mt","mean percent.mt", "min percent.mt", "clusters unfit", 
                                        "num cells after"))

sample_list <- list.dirs(input_path, full.names = F, recursive = F)
for (i in 1:length(sample_list)) { #usually: 3:length(sample_list)
  print(sample_list[i])
  dir.create(paste(output_path, sample_list[i], sep = ""))
  # load sample object
  seurat_obj <- readRDS(file = paste(input_path, sample_list[i], "/", sample_list[i], ".rds", sep = ""))
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
  num_after <- dim(seurat_obj)[2]
  toRemove <- paste(toRemove, collapse=', ' ) # collapse to create one 'string' of cluster numbers
  # check for #cells after removal 
  col2add <- c(col2add, toRemove, num_after)
  # insert relevant data to summary table
  summary_DF[[sample_list[[i]]]] <- col2add
  saveRDS(seurat_obj, file = paste(output_path, sample_list[i],"/", sample_list[i], ".rds", sep = ""))
  rm(j,k,toRemove, temp_df, col2add, seurat_obj,num_after)
}

summary_DF <- apply(summary_DF, 2, as.character)
row.names(summary_DF) <- c("num cells before", 
                           "max nFeature", "mean nFeature", "min nFeature", 
                           "max nCounts","mean nCounts", "min nCounts", 
                           "max percent.mt","mean percent.mt", "min percent.mt","clusters unfit",
                           "num cells after")
write.csv(summary_DF, file = paste(output_path, "summary_DF.csv", sep = ""), row.names = TRUE)

#send into basic analysis again