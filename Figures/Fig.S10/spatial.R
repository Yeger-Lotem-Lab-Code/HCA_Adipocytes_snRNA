

library(Seurat)
library(ggplot2)

##Download the visium rds object and the section tissue position data from Mendeley data at: https://data.mendeley.com/datasets/3bs5f8mvbs/1
seurat_object <- readRDS("/Users/mayaziv/Downloads/se-object.visium_baseline.rds")
tissue_positions <- read.table("/Users/mayaziv/Downloads/tissue_positions_list_s44.csv", header = FALSE, sep = ",")

colnames(tissue_positions) <- c("barcode","v2", "array_row", "array_col", "pos_x", "pos_y")
cells_to_keep <- WhichCells(seurat_object, expression = novaseq_id == "S44")
seurat_subset <- subset(seurat_object, cells = cells_to_keep)
Idents(seurat_subset) <- seurat_subset@meta.data$cluster_anno

seurat_barcodes <- colnames(seurat_subset)

colnames(seurat_subset) <- gsub("_2$", "", colnames(seurat_subset))
seurat_subset@meta.data$barcode <- rownames(seurat_subset@meta.data)
spatial_metadata <- merge(seurat_subset@meta.data, tissue_positions,  by.x = "barcode", by.y = "barcode", all.x = TRUE)
head(spatial_metadata)


# Assuming spatial_metadata is already merged as shown in previous steps
seurat_subset[["pos_x"]] <- spatial_metadata$pos_x
seurat_subset[["pos_y"]] <- spatial_metadata$pos_y

max_x <- max(seurat_subset[["pos_x"]], na.rm = TRUE)
max_y <- max(seurat_subset[["pos_y"]], na.rm = TRUE)
seurat_subset[["pos_x_scaled"]] <- seurat_subset[["pos_x"]] / max_x
seurat_subset[["pos_y_scaled"]] <- seurat_subset[["pos_y"]] / max_y

Adipocyte_subset <- subset(seurat_subset, idents =c("Adipocyte 1", "Adipocyte 2", "Adipocyte 3"), invert= FALSE)
Pre_Adipocyte_subset <- subset(seurat_subset, idents =c("Adipocyte 1", "Adipocyte 2", "Adipocyte 3", "Preadipocyte 1", "Preadipocyte 2", "Preadipocyte 3", "Preadipocyte 4"), invert= FALSE)


# Extracting metadata and PTPRC expression data
data_to_plot <- FetchData(seurat_subset, vars = c("pos_x_scaled", "pos_y_scaled", "PTPRC", "ANK2", "PTPRB", "PDE4D",  "PLIN1", "CD36", "cluster_anno"))
adipo_data <- FetchData(Adipocyte_subset, vars = c("pos_x_scaled", "pos_y_scaled",  "PTPRC", "ANK2", "PTPRB", "PDE4D",  "PLIN1", "CD36", "cluster_anno"))
pre_adipo_data <- FetchData(Pre_Adipocyte_subset, vars = c("pos_x_scaled", "pos_y_scaled",  "PTPRC", "ANK2", "PTPRB", "PDE4D",  "PLIN1", "CD36", "cluster_anno"))

markers <- c("PTPRC",  "ANK2", "PTPRB", "PDE4D", "PLIN1", "CD36")  # Add more markers if needed
sample="S44"
for (marker in markers) {
  print(marker)
  marker= "PTPRC"
  # custom_colors <- c("Adipocyte 1" = "#dfc27d", "Adipocyte 2" = "#80cdc1", "Adipocyte 3" = "#018571")
  
  tissue_plot <- ggplot(adipo_data, aes(x = pos_x_scaled, y = pos_y_scaled, color = PTPRC, shape = `cluster_anno`)) +
    geom_point(size = 2, stroke = 1.5) +  # Plot points
    scale_color_gradient(low = "#fee5d9", high = "#cb181d") +  # Color gradient for gene expression
    scale_shape_manual(values = 1:20) +
    # Define a more extensive set of shapes if you have many clusters
    #scale_shape_manual(values = custom_shapes)
    theme_minimal() +  # Minimal theme
    labs(title = paste("Cell Distribution by " ,marker, " Expression"),
         x = "X Coordinate",
         y = "Y Coordinate",
         color = paste(marker," expression"),
         shape = "Cluster annotation") +
    theme(legend.position = "right")
  jpeg(filename = paste(output_path, "tissue_plot_adipo", marker,".jpeg", sep = ""),
       width = 1500,  # Width in pixels
       height = 1500,  # Height in pixels
       res = 150)
  print(tissue_plot)
  dev.off()  
  
  }