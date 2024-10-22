suppressMessages(library(Seurat))
suppressMessages(library(CellChat))
suppressMessages(library(patchwork))
options(stringsAsFactors = FALSE)
suppressMessages(library(future))
suppressMessages(library(ggalluvial))
suppressMessages(library(NMF))

#upload your seurat object
seurat_object <- readRDS(paste0(input_path, "adipo10.rds"))

# # Extract the CellChat input files from a Seurat V3 object

data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(seurat_object)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
# Create a CellChat object
cellchat_vis11 <- createCellChat(object = data.input, meta = meta, group.by = "group")
rm(seurat_object, data.input, labels, meta)

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# showDatabaseCategory(CellChatDB)
# dplyr::glimpse(CellChatDB$interaction)
# choose whether to use all of cellchatDB or just a subset:
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- CellChatDB #use all cellchatDB
cellchat_vis11@DB <- CellChatDB.use
# 
# # Preprocessing the expression data for cell-cell communication analysis
cellchat_vis11 <- subsetData(cellchat_vis11) # This step is necessary even if using the whole database

# 
cellchat_vis11 <- identifyOverExpressedGenes(cellchat_vis11)
cellchat_vis11 <- identifyOverExpressedInteractions(cellchat_vis11)
# # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in
# # the function `computeCommunProb()` in order to use the projected data)
 cellchat_vis11 <- projectData(cellchat_vis11, PPI.human)
# 
# # Compute the communication probability and infer cellular communication network
cellchat_vis11 <- computeCommunProb(cellchat_vis11, population.size = TRUE)
# # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_vis11 <- filterCommunication(cellchat_vis11, min.cells = 10)
# 
# #Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat_vis11, slot.name= "netP")
write.csv(df.net, file = paste0(output_path, "df.adipo_sub_path_cellchat_adipo_vis_5_8_netp.csv"))
df.net <- subsetCommunication(cellchat_vis11)
# 
saveRDS(cellchat_vis11, file = paste0(output_path, "cellchat_adipo_vis_5_8.rds"))
write.csv(df.net, file = paste0(output_path, "df.adipo_sub_path_cellchat_adipo_vis_5_8.csv"))
rm(CellChatDB, df.net)


input_path = "insert_your_input_path"
output_path = "insert_your_output_path"

cellchat_adipo <- readRDS(file = paste0(input_path, "cellchat_adipo_vis_5_8.rds"))
df.net_adipo <- read.csv(file = paste0(input_path, "df.adipo_sub_path_cellchat_adipo_vis_5_8.csv"), header = T, row.names = 1)

cellchat_adipo <- computeCommunProbPathway(cellchat_adipo)
labels.levels= c("VA1","VA2","VA3", "VA4", "VA5", "VA6", "VA7", "VA8", "V_Mac1", "V_Mac2", "V_Monocytes","V_Pre-Myeloid/Mac","V_Pre-Myeloid","V_S_Mast cells", "V_T cells/NK", "V_B cells")

cellchat_adipo <-updateClusterLabels(
  cellchat_adipo,
  old.cluster.name = NULL,
  new.cluster.name = NULL,
  new.order = labels.levels,
  new.cluster.metaname = "new.labels"
)
#levels(cellchat_adipo@meta$labels)[levels(cellchat_adipo@meta$labels) =='V_Monocytes' ] <- 'Monocytes' 



a <- 'ADIPONECTIN'


cellchat_adipo <- aggregateNet(cellchat_adipo)
a="ADIPONECTIN"

#visualize:

## FIGURE 7:
# Compute the network centrality scores
cellchat_adipo <- netAnalysis_computeCentrality(cellchat_adipo, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups

pdf(file = paste(output_path , "Aadipoq_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'ADIPONECTIN', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()

pdf(file = paste(output_path ,a, "netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = a, width = 8,
                                  height = 2.5, font.size = 10)
dev.off()

pdf(file = paste(output_path , "CD45_netAnalysis_signalingRole_network.pdf"))



netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'CD45', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()

pdf(file = paste(output_path , "cd22_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'CD22', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()

pdf(file = paste(output_path , "laminin_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'LAMININ', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()
pdf(file = paste(output_path , "Adgre5_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'ADGRE5', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()
pdf(file = paste(output_path , "Angptl_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'ANGPTL', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()
pdf(file = paste(output_path , "Annexin_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'ANNEXIN', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()
pdf(file = paste(output_path , "vegf_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'VEGF', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()
pdf(file = paste(output_path , "ANGPT_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'ANGPT', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()
pdf(file = paste(output_path , "collagen_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'COLLAGEN', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()

pdf(file = paste(output_path , "VISFATIN_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'VISFATIN', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()

pdf(file = paste(output_path , "SEMA3_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'SEMA3', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()
pdf(file = paste(output_path , "leptin_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'LEP', width = 8,
                                  height = 2.5, font.size = 10)
dev.off()
pdf(file = paste(output_path , "il16_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'IL16', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()


pdf(file = paste(output_path , "Annexin_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'ANNEXIN', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()
pdf(file = paste(output_path , "COLLAGEN_netAnalysis_signalingRole_network.pdf"))

#jpeg(paste(output_path ,a, "netAnalysis_signalingRole_network.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'COLLAGEN', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()
pdf(file = paste(output_path , "sema3_netAnalysis_signalingRole_network.pdf"))


netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'SEMA3', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()


##Until here FIG.7

## FIG S21-ANNEXIN signaling network in hSAT and hVAT:


#visualize seperatly:
 mat <- cellchat_adipo@net$weight
# par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  pdf(file = paste(output_path , i,".pdf",  sep = ""))
  jpeg(paste(output_path , i,".jpeg", sep = ""),
  width = 1000, height = 750)
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.size = 0.05,
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

# # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')

pdf(file = paste(output_path ,"mastcells_heatmap.pdf"))
netVisual_heatmap(cellchat_adipo,sources.use="Mast cells", color.heatmap = "Reds")
dev.off()
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
pdf(file = paste(output_path ,"mac2_adipocytres_chord signaling cells.pdf"))
netVisual_chord_gene(cellchat_adipo, targets.use = c("SA1","SA2","SA3", "SA4", "SA5","SA6","SA7"), sources.use = "S_Mac2", lab.cex=1)
dev.off()

jpeg(paste(output_path , "netvisual.jpeg", sep = ""),
     width = 1000, height = 750)
netVisual_barplot(cellchat_adipo,comparison=c("VA1", "VA2"), title.name = "netVisual_barplot")
dev.off()
print (cellchat_adipo@netP$pathways)
par(mfrow = c(1,1))



jpeg(paste(output_path , a,"_vertex.receiver.jpeg", sep = ""),
 width = 1000, height = 750)
netVisual_aggregate(cellchat_adipo, signaling = a,  vertex.receiver = seq(1,4))
dev.off()
# 
pdf(file = paste(output_path ,a,"_plotGeneExpression.pdf"))
plotGeneExpression(cellchat_adipo, signaling = a)
dev.off()
# 
pdf(file = paste(output_path ,a,"_circle.pdf"))
netVisual_aggregate(cellchat_adipo, signaling = a, layout = "circle")
dev.off()

#adjust colors to your labels
# cols_mapped_adipo <- c(
#   "Adipocytes 1" = "#6600CC",
#   "Adipocytes 2" = "#6600CC",
#   "Adipocytes 3" = "#6600CC",
#   "Macrophages 1" = "#CC0066",
#   "Macrophages 2" = "#CC0066",
#   "Macrophages 3" = "#CC0066",
#   "Monocytes" = "#ff7f00",
#   "S_Mast cells" = "#FFB266",
#   "Pre-Myeloid" = "#BC8F8F",
#   "B cells" = "#CD853F",
#   "T.NK" = "#00CC00"
# )

pdf(file = paste(output_path ,"cd45_mast_source_color_chord.pdf"))
netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
dev.off()
pdf(file = paste(output_path ,"laminin_mast_source_color_chord.pdf"))
#jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
#     width = 1000, height = 750)
#netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
netVisual_aggregate(cellchat_adipo,sources.use="V_S_Mast cells", signaling = 'LAMININ', layout = "chord")
dev.off()
# pdf(file = paste(output_path ,"annexin_mast_source_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo,sources.use="V_S_Mast cells", signaling = 'ANNEXIN', layout = "chord")
# dev.off()
# pdf(file = paste(output_path ,"Visfatin_mast_source_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo,sources.use="V_S_Mast cells", signaling = 'VISFATIN', layout = "chord")
# dev.off()
# 
# pdf(file = paste(output_path ,"sema3_mast_source_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo,sources.use="V_S_Mast cells", signaling = 'SEMA3', layout = "chord")
# dev.off()
# pdf(file = paste(output_path ,"collagen_mast_source_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo,sources.use="V_S_Mast cells", signaling = 'COLLAGEN', layout = "chord")
# dev.off()
# 
# pdf(file = paste(output_path ,"angptl_mast_source_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo,sources.use="V_S_Mast cells", signaling = 'ANGPTL', layout = "chord")
# dev.off()
# 
# pdf(file = paste(output_path ,"CD45_mast_source_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo,sources.use="V_S_Mast cells", signaling = 'CD45', layout = "chord")
# dev.off()
# pdf(file = paste(output_path ,"visfatin_heatmap.pdf"))
# 
# #jpeg(paste(output_path ,a, "_heatmap.jpeg", sep = ""),
# #    width = 1000, height = 750)
# netVisual_heatmap(cellchat_adipo, signaling = 'VISFATIN', color.heatmap = "Reds")
# dev.off()
# pdf(file = paste(output_path ,"lep_heatmap.pdf"))
# 
# #jpeg(paste(output_path ,a, "_heatmap.jpeg", sep = ""),
# #    width = 1000, height = 750)
# netVisual_heatmap(cellchat_adipo, signaling = 'LEP', color.heatmap = "Reds")
# dev.off()
# pdf(file = paste(output_path ,"SEMA3_heatmap.pdf"))
# 
# #jpeg(paste(output_path ,a, "_heatmap.jpeg", sep = ""),
# #    width = 1000, height = 750)
# netVisual_heatmap(cellchat_adipo, signaling = 'SEMA3', color.heatmap = "Reds")
# dev.off()
# pdf(file = paste(output_path ,"annexin_heatmap.pdf"))
# 
# #jpeg(paste(output_path ,a, "_heatmap.jpeg", sep = ""),
# #    width = 1000, height = 750)
# netVisual_heatmap(cellchat_adipo, signaling = 'ANNEXIN', color.heatmap = "Reds")
# dev.off()

# pdf(file = paste(output_path ,"il16_heatmap.pdf"))
# 
# #jpeg(paste(output_path ,a, "_heatmap.jpeg", sep = ""),
# #    width = 1000, height = 750)
# netVisual_heatmap(cellchat_adipo, signaling = 'IL16', color.heatmap = "Reds")
# dev.off()
# pdf(file = paste(output_path ,"leptin_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo, signaling = 'LEP', layout = "chord")
# # dev.off()
# pdf(file = paste(output_path ,"annexin_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo, signaling = 'ANNEXIN', layout = "chord")
# dev.off()

# pdf(file = paste(output_path ,"leptin_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo, signaling = 'LEP', layout = "chord")
# dev.off()

pdf(file = paste(output_path ,"adipoq_color_chord.pdf"))
#jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
#     width = 1000, height = 750)
#netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
netVisual_aggregate(cellchat_adipo, signaling = 'ADIPONECTIN', layout = "chord")
dev.off()

# pdf(file = paste(output_path ,"IL16_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo, signaling = 'IL16', layout = "chord")
# dev.off()

#pdf(file = paste(output_path ,"VEGF_color_chord.pdf"))
#jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
#     width = 1000, height = 750)
#p1 <- netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'VEGF', layout = "chord", pt.title = 15,title.space	=1)
jpeg(paste0(input_path, "chord_mast_vegf.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'VEGF', layout = "chord", pt.title = 15,title.space	=1)
dev.off()

jpeg(paste0(input_path, "chord_mast_laminin.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'LAMININ', layout = "chord", pt.title = 15,title.space	=1)
dev.off()

jpeg(paste0(input_path, "chord_mast_angpt.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'ANGPT', layout = "chord", pt.title = 15,title.space	=1)
dev.off()

jpeg(paste0(input_path, "chord_mast_collagen.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'COLLAGEN', layout = "chord", pt.title = 15,title.space	=1)
dev.off()

jpeg(paste0(input_path, "chord_mast_sema3.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'SEMA3', layout = "chord", pt.title = 15,title.space	=1)
dev.off()

jpeg(paste0(input_path, "chord_mast_ANNEXIN.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'ANNEXIN', layout = "chord", pt.title = 15,title.space	=1)
dev.off()

jpeg(paste0(input_path, "chord_mast_IL16.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'IL16', layout = "chord", pt.title = 15,title.space	=1)
dev.off()

jpeg(paste0(input_path, "chord_mast_vegf.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'VEGF', layout = "chord", pt.title = 15,title.space	=1)
dev.off()


p1 <- netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'LAMININ', layout = "chord", pt.title = 15,title.space	=1)
jpeg(paste0(input_path, "chord_mast_LAMININ.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
print(p1)
dev.off()

p1 <- netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'ANGPT', layout = "chord", pt.title = 15,title.space	=1)
jpeg(paste0(input_path, "chord_mast_ANGPT.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
print(p1)
dev.off()

p1 <- netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'COLLAGEN', layout = "chord", pt.title = 15,title.space	=1)
jpeg(paste0(input_path, "chord_mast_collagen.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
print(p1)
dev.off()

p1 <- netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'SEMA3', layout = "chord", pt.title = 15,title.space	=1)
jpeg(paste0(input_path, "chord_mast_SEMA3.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
print(p1)
dev.off()

p1 <- netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'IL16', layout = "chord", pt.title = 15,title.space	=1)
jpeg(paste0(input_path, "chord_mast_il16.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
print(p1)
dev.off()

p1 <- netVisual_aggregate(cellchat_adipo, sources.use="S_Mast cells", signaling = 'ANNEXIN', layout = "chord", pt.title = 15,title.space	=1)
jpeg(paste0(input_path, "chord_mast_ANNEXIN.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
print(p1)
dev.off()



# pdf(file = paste(output_path ,"il16_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo,targets.use="SA3", signaling = 'IL16', layout = "chord")
# dev.off()

# pdf(file = paste(output_path ,"angpt_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo,targets.use="SA2", signaling = 'ANGPT', layout = "chord")
# dev.off()
# 
# pdf(file = paste(output_path ,"vegf_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo,targets.use="SA2", signaling = 'VEGF', layout = "chord")
# dev.off()
# pdf(file = paste(output_path ,"laminin_color_chord.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# #netVisual_aggregate(cellchat_adipo,color.use= cols_mapped_adipo,  sources.use="S_Mast cells"	, signaling = a, layout = "chord")
# netVisual_aggregate(cellchat_adipo,targets.use="SA5", signaling = 'LAMININ', layout = "chord")
# dev.off()


# pdf(file = paste(output_path ,"va1_all_chord between cells.pdf"))
# #netVisual_chord_gene(cellchat_adipo, source.use = c('T.NK'), targets.use = c("VA3"), lab.cex = 0.5,legend.pos.y = 30)
# netVisual_chord_gene(cellchat_adipo, sources.use = 'VA1', targets.use = c('VA2','VA3','VA4','VA5','VA6','VA7','VA8'), lab.cex = 0.5,legend.pos.y = 30)
# dev.off()
# pdf(file = paste(output_path ,"t2_va3_chord between cells.pdf"))
# #netVisual_chord_gene(cellchat_adipo, source.use = c('T.NK'), targets.use = c("VA3"), lab.cex = 0.5,legend.pos.y = 30)
# netVisual_chord_gene(cellchat_adipo, sources.use = 'S_T cells 2', targets.use = c('SA3','SA4'), lab.cex = 0.5,legend.pos.y = 30)
# dev.off()
# pdf(file = paste(output_path ,"mast_va2_chord between cells.pdf"))
# netVisual_chord_gene(cellchat_adipo, source.use = c('V_S_Mast cells'), targets.use = c("VA2"), lab.cex = 0.5,legend.pos.y = 30)
# dev.off()
# pdf(file = paste(output_path ,"mast_va3_chord between cells.pdf"))
# netVisual_chord_gene(cellchat_adipo, source.use = c('V_S_Mast cells'), targets.use = c("VA3"), lab.cex = 0.5,legend.pos.y = 30)
# dev.off()
# pdf(file = paste(output_path ,"va1_all_va3_chord between cells.pdf"))
# netVisual_chord_gene(cellchat_adipo, source.use = c('VA1'), targets.use = c("VA2","VA3","VA4","VA5","VA6", "VA7", "VA8"), lab.cex = 0.5,legend.pos.y = 30)
# dev.off()

# pdf(file = paste(output_path ,"color_chord_sa2_source.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# netVisual_aggregate(cellchat_adipo, sources.use="SA2"	, signaling = c('VEGF'), layout = "chord")
# dev.off()
# 
# pdf(file = paste(output_path ,"color_chord_sa2_target.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# netVisual_aggregate(cellchat_adipo  ,pt.title = 3,
#                     title.space = 1,targets.use="SA2"	, signaling = c('VEGF'), layout = "chord")
# dev.off()
# pdf(file = paste(output_path ,"color_chord_sa3_target_annexin_il16.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# netVisual_aggregate(cellchat_adipo,  target.use="SA3"	, signaling = c('IL16'), layout = "chord")
# dev.off()
# pdf(file = paste(output_path ,"color_chord_sa3_source_annexin_il16.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# netVisual_aggregate(cellchat_adipo,  source.use="SA3"	, signaling = c('ANNEXIN'), layout = "chord")
# dev.off()
# pdf(file = paste(output_path ,"color_chord_sa1_source_adipo_lep.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# netVisual_aggregate(cellchat_adipo,  source.use="SA1"	, signaling = c('ADIPONECTIN'), layout = "chord")
# dev.off()
# pdf(file = paste(output_path ,"color_chord_sa1_target_adipo_lep.pdf"))
# #jpeg(paste(output_path , a,"color_chord.jpeg", sep = ""),
# #     width = 1000, height = 750)
# netVisual_aggregate(cellchat_adipo,  target.use="SA1"	, signaling = c('LEP'), layout = "chord")
# dev.off()
# 
# pdf(file = paste(output_path ,a,"netAnalysis_contribution.pdf"))
# 
# #jpeg(paste(output_path , a,"netAnalysis_contribution.jpeg", sep = ""),
# #     width = 1000, height = 750)
# netAnalysis_contribution(cellchat_adipo, signaling = a)
# dev.off()




library(ggalluvial)
library(NMF)
# Identify and visualize outgoing communication pattern of secreting cells
selectK(cellchat_adipo, pattern = "outgoing")
nPatterns = 4
#jpeg(paste(output_path , "identifyCommunicationPatterns_outgoing.jpeg", sep = ""))
pdf(file = paste(output_path ,"identifyCommunicationPatterns_outgoing.pdf"),height=20)
cellchat_adipo <- identifyCommunicationPatterns(cellchat_adipo, pattern = "outgoing", k = nPatterns, height=20)
dev.off()
# river plot
pdf(file = paste(output_path ,"river plot_outgoing.pdf"))

#jpeg(paste(output_path , "river plot_outgoing.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_river(cellchat_adipo, pattern = "outgoing")
dev.off()
# dot plot
pdf(file = paste(output_path ,"dot_plot_outgoing.pdf"))

#jpeg(paste(output_path , "dot_plot_outgoing.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_dot(cellchat_adipo, pattern = "outgoing")
dev.off()

selectK(cellchat_adipo, pattern = "incoming")
nPatterns = 4
# jpeg(paste(output_path , "identifyCommunicationPatterns_incoming.jpeg", sep = ""))
#par(cex=0.5, cex.main=0.5, cex.lab = 0.5, cex.sub=0.8)
pdf(file = paste(output_path ,"identifyCommunicationPatterns_incoming.pdf"), height=20)
cellchat_adipo <- identifyCommunicationPatterns(cellchat_adipo, pattern = "incoming", k = nPatterns, height=20)
dev.off()
# river plot
pdf(file = paste(output_path ,"river_plot_incoming.pdf"))

#jpeg(paste(output_path , "river_plot_incoming.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_river(cellchat_adipo, pattern = "incoming")
dev.off()
# dot plot
pdf(file = paste(output_path ,"dot_plot_incoming.pdf"))

#jpeg(paste(output_path , "dot_plot_incoming.jpeg", sep = ""),
#    width = 1000, height = 750)
netAnalysis_dot(cellchat_adipo, pattern = "incoming")
dev.off()
# Manifold and classification learning analysis of signaling networks
# functional similarity
cellchat_adipo <- computeNetSimilarity(cellchat_adipo, type = "functional")
library(reticulate)
# library(umap) #solve this. needs python package: umap-learn

cellchat_adipo <- netEmbedding(cellchat_adipo, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat_adipo <- netClustering(cellchat_adipo, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
pdf(file = paste(output_path ,"2D_signaling_network_finctional.pdf"))

#jpeg(paste(output_path , "2D_signaling_network_finctional.jpeg", sep = ""),
#    width = 1000, height = 750)
netVisual_embedding(cellchat_adipo, type = "functional", label.size = 3.5)
dev.off()
# structure similarity
cellchat_adipo <- computeNetSimilarity(cellchat_adipo, type = "structural")
cellchat_adipo <- netEmbedding(cellchat_adipo, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat_adipo <- netClustering(cellchat_adipo, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
pdf(file = paste(output_path ,"2D_signaling_network_structural.pdf"))

#jpeg(paste(output_path , "2D_signaling_network_structural.jpeg", sep = ""),
#    width = 1000, height = 750)
netVisual_embedding(cellchat_adipo, type = "structural", label.size = 3.5)
dev.off()
