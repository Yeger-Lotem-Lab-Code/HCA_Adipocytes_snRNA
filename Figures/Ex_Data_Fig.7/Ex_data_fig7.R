## Extended data fig 7

##Addition o the code for Main fig 7 in "cellchat.R"
##Run twice for each immune cells+adipocytes cellchat objects once for hSAt and once for hVAT

jpeg(paste0(input_path, "chord_mast_sema3.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
netVisual_aggregate(cellchat_adipo, sources.use="Mast cells", signaling = 'SEMA3', layout = "chord", pt.title = 15,title.space	=1)
dev.off()
jpeg(paste0(input_path, "chord_mast_collagen.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
netVisual_aggregate(cellchat_adipo, sources.use="Mast cells", signaling = 'COLLAGEN', layout = "chord", pt.title = 15,title.space	=1)
dev.off()

jpeg(paste0(input_path, "chord_mast_laminin.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
netVisual_aggregate(cellchat_adipo, sources.use="Mast cells", signaling = 'LAMININ', layout = "chord", pt.title = 15,title.space	=1)
dev.off()

jpeg(paste0(input_path, "chord_mast_IL16.jpeg"), width = 1500, height = 1500, quality = 100, res=300)
netVisual_aggregate(cellchat_adipo, sources.use="Mast cells", signaling = 'IL16', layout = "chord", pt.title = 15,title.space	=1)
dev.off()
