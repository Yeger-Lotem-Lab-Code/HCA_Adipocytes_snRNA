suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))




output_path = "/gpfs0/estiyl/users/mziv/immune/samples_29.1/raw_objects/"
Bender_data_path = "/gpfs0/estiyl/users/mziv/immune/samples_29.1/cellbender_outputs/"
Ranger_data_path = "/gpfs0/estiyl/users/mziv/"

###---load each sample, from cellBender and cellRanger. intersect and then create seurat obj
create_raw_obj <-  function(sample_name, feature_cutoff){
  dir.create(paste(output_path, sample_name, sep = ""))
  Ranger_data <- Read10X_h5(filename = paste(Ranger_data_path, sample_name, "/outs/", "filtered_feature_bc_matrix.h5", 
                                             sep = ""), use.names = TRUE, unique.features = TRUE)
  # Bender_filtered_barcodes <- read.csv(file = paste(Bender_data_path, sample_name, "/", 
  #                                                   sample_name, "_cell_barcodes.csv", sep = ""), header = F)
  Bender_data <- Read10X_h5(file = paste0(Bender_data_path, sample_name, "/", sample_name, "_filtered.h5"), 
                            use.names = TRUE, unique.features = TRUE)
  cells_inBoth <- intersect(colnames(Ranger_data), colnames(Bender_data))
  obj_seurat = CreateSeuratObject(counts = Bender_data, project = sample_name)
  obj_seurat <- subset(obj_seurat, cells = cells_inBoth)
  obj_seurat[["percent.mt"]] <- PercentageFeatureSet(obj_seurat, pattern = "^MT-")
  
  plot1 <- VlnPlot(obj_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  jpeg(paste(output_path,"/", sample_name,"/", sample_name, "_QCmetrics.jpeg", sep = ""), 
       width = 500, height = 500)
  print(plot1)
  dev.off()
  # plot2 <- FeatureScatter(obj_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")+
  #   theme(legend.position = 'none')
  # plot3 <- FeatureScatter(obj_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  #   theme(legend.position = 'none')
  # jpeg(paste(output_path,"/", sample_name,"/", sample_name, "_FeatureScatter.jpeg", sep = ""), 
  #      width = 1000, height = 500)
  # print(plot2 + plot3)
  # dev.off()
  rm(plot1,Bender_data, Ranger_data, cells_inBoth)
  
  #obj_seurat <- subset(obj_seurat, subset = nFeature_RNA > feature_cutoff & percent.mt < 20)
  obj_seurat <- subset(obj_seurat, subset = nFeature_RNA > feature_cutoff )
  saveRDS(obj_seurat, file = (paste(output_path, "/", sample_name,"/", sample_name,".rds", sep="")))
  rm(obj_seurat)
}

sample_list <- list.dirs(Bender_data_path, full.names = F, recursive = F)
for (i in 1:length(sample_list)) 
  {
    print(sample_list[i])
    create_raw_obj(sample_name = sample_list[i], feature_cutoff = 200)
  }
  
