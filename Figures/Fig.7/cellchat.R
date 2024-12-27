suppressMessages(library(Seurat))
suppressMessages(library(CellChat))
suppressMessages(library(patchwork))
options(stringsAsFactors = FALSE)
suppressMessages(library(future))
suppressMessages(library(ggalluvial))
suppressMessages(library(NMF))


 seurat_object <- readRDS(paste0(input_path, "adipo10.rds"))

# # Extract the CellChat input files from a Seurat V3 object

# data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
# labels <- Idents(seurat_object)
# meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
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
# Preprocessing the expression data for cell-cell communication analysis
cellchat_vis11 <- subsetData(cellchat_vis11) # This step is necessary even if using the whole database
# future::plan("multiprocess", workers = 4) # do parallel
plan(cluster)
options(future.globals.maxSize = 10000 * 1024^2)
plan()
availableCores()
# 
cellchat_vis11 <- identifyOverExpressedGenes(cellchat_vis11)
cellchat_vis11 <- identifyOverExpressedInteractions(cellchat_vis11)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in
# the function `computeCommunProb()` in order to use the projected data)
cellchat_vis11 <- projectData(cellchat_vis11, PPI.human)
# 
# # Compute the communication probability and infer cellular communication network
cellchat_vis11 <- computeCommunProb(cellchat_vis11, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_vis11 <- filterCommunication(cellchat_vis11, min.cells = 10)

# #Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat_vis11, slot.name= "netP")
write.csv(df.net, file = paste0(output_path, "df.adipo_sub_path_cellchat_adipo_vis_5_8_netp.csv"))
df.net <- subsetCommunication(cellchat_vis11)
# 
saveRDS(cellchat_vis11, file = paste0(output_path, "cellchat_adipo_vis_5_8.rds"))
write.csv(df.net, file = paste0(output_path, "df.adipo_sub_path_cellchat_adipo_vis_5_8.csv"))
rm(CellChatDB, df.net)

input_path = "/gpfs0/estiyl/users/mziv/Adipose_new_w_o_sc3369/sub_vis_merge/cellchat/31_7/vis/immune_adipo//"
output_path = "/gpfs0/estiyl/users/mziv/Adipose_new_w_o_sc3369/sub_vis_merge/cellchat/31_7/vis/immune_adipo/"


cellchat_adipo <- readRDS(file = paste0(input_path, "cellchat_adipo_vis_5_8.rds"))
df.net_adipo <- read.csv(file = paste0(input_path, "df.adipo_sub_path_cellchat_adipo_vis_5_8.csv"), header = T, row.names = 1)
# continue:
cellchat_adipo <- computeCommunProbPathway(cellchat_adipo)
labels.levels= c("VA1","VA2","VA3", "VA4", "VA5", "VA6", "VA7", "VA8", "V_Mac1", "V_Mac2", "V_Monocytes","V_Pre-Myeloid/Mac","V_Pre-Myeloid","V_S_Mast cells", "V_T cells/NK", "V_B cells")

cellchat_adipo <-updateClusterLabels(
  cellchat_adipo,
  old.cluster.name = NULL,
  new.cluster.name = NULL,
  new.order = labels.levels,
  new.cluster.metaname = "new.labels"
)



a <- 'ADIPONECTIN'


cellchat_adipo <- aggregateNet(cellchat_adipo)
a="ADIPONECTIN"

#visualize:


# # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# # show all the interactions sending from Inflam.FIB
#
# #
# #

# show all the significant interactions (L-R pairs) associated with certain signaling pathways

print (cellchat_adipo@netP$pathways)

# pdf(file = paste(output_path ,"annexin_heatmap.pdf"))
# 
# #jpeg(paste(output_path ,a, "_heatmap.jpeg", sep = ""),
# #    width = 1000, height = 750)
# netVisual_heatmap(cellchat_adipo, signaling = 'ANNEXIN', color.heatmap = "Reds")
# dev.off()





# Compute the network centrality scores
cellchat_adipo <- netAnalysis_computeCentrality(cellchat_adipo, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups

pdf(file = paste(output_path , "Adipoq_netAnalysis_signalingRole_network.pdf"))

netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'ADIPONECTIN', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()

pdf(file = paste(output_path , "LEPTIN_netAnalysis_signalingRole_network.pdf"))
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'LEP', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()


pdf(file = paste(output_path , "ANGPT_netAnalysis_signalingRole_network.pdf"))
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'ANGPT', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()

pdf(file = paste(output_path , "VEGF_netAnalysis_signalingRole_network.pdf"))
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'VEGF', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()

pdf(file = paste(output_path , "CD45_netAnalysis_signalingRole_network.pdf"))
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'CD45', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()

pdf(file = paste(output_path , "Annexin_netAnalysis_signalingRole_network.pdf"))
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'ANNEXIN', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()

pdf(file = paste(output_path , "IL16_netAnalysis_signalingRole_network.pdf"))
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'IL16', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()

pdf(file = paste(output_path , "LAMININ_netAnalysis_signalingRole_network.pdf"))
netAnalysis_signalingRole_network(cellchat_adipo, signaling = 'LAMININ', width = 14,
                                  height = 2.5, font.size = 10)
dev.off()

