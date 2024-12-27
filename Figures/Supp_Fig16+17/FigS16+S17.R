#### fig S16:----
# Sequencing saturation plots per sample from cellRanger

#### fig S17 sampling coverage:----

# figure S17A: hSAT whole tissue:
jpeg(paste0(output_path, "boxplot_sc_5_All_clusters.jpg"), width = 800, height = 600)
boxplot(hVSAT_result,
        main = paste0("Boxplot for all clusters"),
        ylab = "Max Fraction",
        outline = FALSE,
        names = c("2879","6379","9879","13379", "16879", "20379", "23879","27379","30879","34379"),
        las = 2,
        boxwex = 0.3,
        ylim = c(0, max(hVSAT_result, na.rm = TRUE)))
dev.off()

# figure S17B: hSAT adipocytes:
jpeg(paste0(output_path, "boxplot_adipo5_All_clusters.jpg"), width = 800, height = 600)
boxplot(adipo5_result,
        main = paste0("Boxplot for all clusters"),
        ylab = "Max Fraction",
        outline = FALSE,
        names = c("205","1205","2205","3205", "4205", "5205", "7205","8205","9205","10205", "11205"),
        las = 2,
        boxwex = 0.3,
        ylim = c(0, max(adipo5_result, na.rm = TRUE)))
dev.off()

# figure S17C: hVAT whole tissue:
jpeg(paste0(output_path, "boxplot_vis_10_All_clusters.jpg"), width = 800, height = 600)
boxplot(hVAT_result,
        main = paste0("Boxplot for all clusters"),
        ylab = "Max Fraction",
        outline = FALSE,
        names = c("13731","20731","27731","34371","41731", "48731", "55731", "62731","69731","76731"),
        las = 2,
        boxwex = 0.3,
        ylim = c(0, max(hVAT_result, na.rm = TRUE)))
dev.off()

# figure S17D: hVAT adipocytes:
jpeg(paste0(output_path, "boxplot_adipo10_All_clusters_Log.jpg"), width = 800, height = 600)
boxplot(adipo10_result,
        main = paste0("Boxplot for all clusters"),
        ylab = "Max Fraction",
        outline = FALSE,
        names = c("460", "1460", "2460", "3460", "4460", "5460", "6460", "7460", "8460", "9460", "10460", "11460", "12460", "13460", "14460"),
        las = 2,
        boxwex = 0.3,
        ylim = c(0, max(adipo10_result, na.rm = TRUE)))
dev.off()


