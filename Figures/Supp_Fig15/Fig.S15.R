##continue of the codde of Fig7
#Run twice once for each tissue, with the immune+adipocytes cellchat objects:
pdf(file = paste(output_path ,"annexin_heatmap.pdf"))

netVisual_heatmap(cellchat_adipo, signaling = 'ANNEXIN', color.heatmap = "Reds")
dev.off()
