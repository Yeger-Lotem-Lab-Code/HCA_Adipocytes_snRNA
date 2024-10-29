library(riverplot)
library(Seurat)
library(rliger)
library(dplyr)

output_path = ".."
input_path= "..."
seurat_obj2 <- readRDS('Sat_adipocytes.rds')
seurat_obj1 <- readRDS('...emont_data/SC/SAT_emont.rds')
#seurat_obj1 <- subset(seurat_obj1,  idents ="VA6", invert= TRUE)
print(table(seurat_obj1@active.ident))
seurat_obj1$ann <- seurat_obj1@active.ident
seurat_obj2$ann <- seurat_obj2@active.ident
print(table(seurat_obj1@active.ident))

# Find common features and perform anchor transfer
common_features <- intersect(rownames(seurat_obj1), rownames(seurat_obj2))
transfer_anchors <- FindTransferAnchors(reference = seurat_obj1, query = seurat_obj2,reduction="rpca",  dims = 1:30, features = common_features)
predictions <- TransferData(anchorset = transfer_anchors, refdata = seurat_obj1$ann, dims = 1:30)
seurat_obj2 <- AddMetaData(seurat_obj2, metadata = predictions$predicted.id, col.name = "predicted_cluster")
write.csv(predictions, file = paste(output_path, "pred_ann_ref_house_sat.csv", sep = ""), row.names = TRUE)

# Create contingency table
contingency_table <- table(seurat_obj2$ann, seurat_obj2$predicted_cluster)
write.csv(contingency_table, file = paste(output_path, "contingency_table_ref_house_all.csv", sep = ""), row.names = TRUE)
# Convert contingency table to data frame
river_data <- as.data.frame(contingency_table) %>%
  rename(source = Var1, target = Var2, freq = Freq)

# Write out river data
write.csv(river_data, file = paste(output_path, "river_data_ref_house.csv", sep = ""), row.names = TRUE)

# Define colors for source clusters (as you already have)
source_clusters <- unique(river_data$source)
# cluster_colors <- c("VA1" = "#d73027",   # Example: Red
#                     "VA2" = "#f46d43",   # Example: Orange
#                     "VA3" = "#fdae61",   # Example: Yellow
#                     "VA4" = "#fee090",   # Example: Light Yellow
#                     "VA5" = "#e0f3f8",   # Example: Light Blue
#                     "VA6" = "#abd9e9",   # Example: Blue
#                     "VA7" = "#74add1",   # Example: Dark Blue
#                     "VA8" = "#4575b4")   # Example: Darker Blue
cluster_colors <- c("SA1" = "#d73027",   # Example: Red
                    "SA2" = "#f46d43",   # Example: Orange
                    "SA3" = "#fdae61",   # Example: Yellow
                    "SA4" = "#fee090",   # Example: Light Yellow
                    "SA5" = "#e0f3f8",   # Example: Light Blue
                    "SA6" = "#abd9e9",   # Example: Blue
                    "SA7" = "#74add1")   
# cluster_colors <- c("hAd1" = "#d73027",   # Example: Red
#                     "hAd2" = "#f46d43",   # Example: Orange
#                     "hAd3" = "#fdae61",   # Example: Yellow
#                     "hAd4" = "#fee090",   # Example: Light Yellow
#                     "hAd5" = "#e0f3f8",   # Example: Light Blue
#                     "hAd6" = "#abd9e9",   # Example: Blue
#                     "hAd7" = "#74add1"  # Example: Dark Blue
#                     )
# Filter cluster colors to only those that exist in the source
cluster_colors <- cluster_colors[source_clusters]

# Create a neutral color for the target nodes (for example, grey)
neutral_color <- "#cccccc"  # You can choose any neutral color you prefer

# Create nodes with sorted order for both source and target
sorted_source <- sort(unique(river_data$source))
sorted_target <- sort(unique(river_data$target))

# nodes <- data.frame(
#   ID = c(sorted_source, sorted_target),
#   x = c(rep(1, length(sorted_source)), rep(2, length(sorted_target))),
#   labels = c(sorted_source, sorted_target),
#   col = c(cluster_colors[match(sorted_source, source_clusters)], rep(neutral_color, length(sorted_target))),
#   stringsAsFactors = FALSE
# )

nodes <- data.frame(
  ID = c(sorted_source, sorted_target),
  x = c(rep(1, length(sorted_source)), rep(2, length(sorted_target))),
  labels = NA,
  col = c(cluster_colors[match(sorted_source, source_clusters)], rep(neutral_color, length(sorted_target))),
  stringsAsFactors = FALSE
)

# Prepare the list of styles for source nodes and apply neutral color for target nodes
styles <- list()

# Define styles for source clusters (with specific colors)
for (cluster in sorted_source) {
  styles[[cluster]] <- list(col = cluster_colors[cluster], lty = 0, textcol = "black", cex.lab = 20, cex = 50, textpos =1)
}
# for (cluster in sorted_target) {
#   styles[[cluster]] <- list( lty = 0, textcol = "black", cex.lab = 20, cex = 50, textpos =2)
# }


# Create edges
edges <- river_data %>%
  rename(N1 = source, N2 = target, Value = freq) 

# Generate the riverplot
#river <- makeRiver(nodes = nodes, edges = edges,  node_labels="")
river <- makeRiver(nodes = nodes, edges = edges, node_styles = styles)

# Save plot to PDF
pdf(paste(output_path, "river_plot_sat_house_emont_ref_10.pdf"), width = 25, height = 40)
plot(river, srt = 90, cex = 1)
dev.off()
