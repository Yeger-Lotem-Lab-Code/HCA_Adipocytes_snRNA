##Fig3_ GSEA_color_based_on python code _"gsea_color.py"


input_path <- "insert_your_input_path"

obj_seurat <- readRDS(paste0(input_path,'adipo10_diet.rds', sep=""))
obj_seurat$adipo_ann <- obj_seurat@active.ident

Idents(obj_seurat) <- obj_seurat$adipo_ann
##read the adjusted gsea score output from the "gsea_color.py"
metadata_df <- read.csv(paste0(output_path,'hvat_adj_score.csv'), header= TRUE)

##add the score of each enrichment pathway to each cells 
for (i in 2:ncol(metadata_df)) {
  cluster_name <- colnames(metadata_df)[i]
  metadata_VAlues <- metadata_df[, i]
  
  for (j in 1:nrow(metadata_df)) {
    cluster <- metadata_df$Cluster[j]
    print(cluster)
    
    # Find the indices of cells with matching annotation
    matching_indices <- which(obj_seurat$adipo_ann == cluster)
    # Update the Seurat object's metadata using the matching indices
    obj_seurat@meta.data[matching_indices, cluster_name] <- metadata_VAlues[j]
  }
}

a="vis"
p1 <- FeaturePlot(obj_seurat, features='Adaptive.immune')+ scale_colour_gradientn( colors= c("#d9d9d9","#ffeda0", "#bd0026"), limits = c(-1.5,1.5 ))
jpeg(paste0(input_path, a,"_adaptive.jpeg"), width = 1700, height = 1500, quality = 100, res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(obj_seurat, features='Angiogenesis')+ scale_colour_gradientn( colors= c("#d9d9d9","#ffeda0", "#bd0026"), limits = c(-1.5,1.5 ))
jpeg(paste0(input_path, a, "_Angiogenesis.jpeg"), width = 1700, height = 1500, quality = 100, res=300)
print(p1)
dev.off()


p1 <- FeaturePlot(obj_seurat, features='Innate.immune')+ scale_colour_gradientn( colors= c("#d9d9d9","#ffeda0", "#bd0026"), limits = c(-1.5,1.5 ))
jpeg(paste0(input_path, a,"_innate.jpeg"), width = 1700, height = 1500, quality = 100, res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(obj_seurat, features='ECM')+ scale_colour_gradientn( colors= c("#d9d9d9","#ffeda0", "#bd0026"), limits = c(-1.5,1.5 ))
jpeg(paste0(input_path, a,"_ecm.jpeg"), width = 1700, height = 1500, quality = 100, res=300)
print(p1)
dev.off()

p1 <- FeaturePlot(obj_seurat, features='Ribosomal')+ scale_colour_gradientn( colors= c("#d9d9d9","#ffeda0", "#bd0026"), limits = c(-1.5,1.5 ))
jpeg(paste0(input_path,a, "_ER_ribosomal.jpeg"), width = 1700, height = 1500, quality = 100, res=300)
print(p1)
dev.off()



p1 <- FeaturePlot(obj_seurat, features='Lipid.Metabolism')+ scale_colour_gradientn( colors= c("#d9d9d9","#ffeda0", "#bd0026"), limits = c(-1.5,1.5 ))
jpeg(paste0(input_path, "hAd_Lipid_metabolism.jpeg"), width = 1700, height = 1500, quality = 100, res=300)
print(p1)
dev.off()
