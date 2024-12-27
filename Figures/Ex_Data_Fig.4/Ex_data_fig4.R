#### EXTENDED FIGURE 4:----
#correlations of only bulk
decon_param <- colnames(cVs_sn_decon[5:61])
cVs_sn_decon$Tissue <- factor(cVs_sn_decon$Tissue)
cVs_sn_decon$cluster <- factor(cVs_sn_decon$cluster)
# dir.create(paste0(output_path, "correlations/"))
corr_table_results <- data.frame(Tissue=NA,param=NA, cluster=NA, rho=NA, p_val=NA)
for (param in decon_param) {     # iterate over param 
  if (typeof(cVs_sn_decon$param)!="character"){
    fig_path <- paste0(output_path, "correlations/", param, "_corr/")
    dir.create(fig_path)   # create sub folder for param corr
    for (tis in levels(cVs_sn_decon$Tissue)) {      # iterate over tis, cluster
      print(tis)
      for (clus in levels(cVs_sn_decon$cluster)) {
        print(clus)
        X <- cVs_sn_decon %>% 
          subset(Tissue==tis & cluster==clus)
        try({
          res1 <- cor.test(X$value,X[,param] ,  method = "spearman",exact=FALSE)
          line2add <- c(tis, param, clus, res1$estimate, res1$p.value)
          corr_table_results <- rbind(corr_table_results, line2add)
          p1 <- ggscatter(X, x = "value", y = param, size=4,
                          add = "reg.line", conf.int = TRUE,
                          cor.coef = F, cor.method = "spearman",
                          xlab = "Percentage of cells", ylab = param)+ 
            geom_smooth(method=lm,fill = "#d9d9d9", color="black")+
            theme(axis.text=element_text(size=12))
          jpeg(paste0(fig_path, tis,"_",clus,"_scatterplot.jpeg"), width = 500, height = 500)
          print(p1)
          dev.off()
        }) # trycatch ends 
      } # iter over cluster
    } # iter over tis
  }} # iter over param ends
# adjust p-val:
corr_table_results$param <- factor(corr_table_results$param)
corr_table_results$Tissue <- factor(corr_table_results$Tissue)
cVs_corr_results2 <- data.frame()
for (tis in levels(corr_table_results$Tissue)){
  for (par in levels(corr_table_results$param)){
    print(par)
    x <- corr_table_results %>% subset(Tissue==tis & param==par) %>% adjust_pvalue(p.col="p_val", method = "BH")
    cVs_corr_results2 <- rbind(cVs_corr_results2, x)
  }
}
cVs_corr_results2 <- cVs_corr_results2[complete.cases(cVs_corr_results2), ]
