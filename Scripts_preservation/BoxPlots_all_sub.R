library(ggplot2)
library(reshape2)

table <-  read.csv(file = paste0("/gpfs0/estiyl/users/tomul/sc_5/1/1.1/maxFractions.csv"))
table <- data.frame(table$maxFractions)
colnames(table) <-c( "table_1.1")
for (i in 1:10) {
  for (j in 1:10) {
    if ((i == 1 && j == 1) || j == 6) {
      next
    }
    tableToAdd <-  read.csv(file = paste0("/gpfs0/estiyl/users/tomul/sc_5/", j, "/", j, ".", i, "/maxFractions.csv"))
    tableToAdd <- data.frame(tableToAdd$maxFractions)
    col_name <- paste0("table_", j, ".", i)
    
    # Add the current table with descriptive column names
    tableTmp <- cbind(table, setNames(tableToAdd, col_name))
    table <- tableTmp
    
  }
}
write.csv(table , paste0("/gpfs0/estiyl/users/tomul/sc_5/results/combinedTable_sc_5.csv"))
# making boxplots for each cluster
  object <-  readRDS(file = paste0("/gpfs0/estiyl/users/tomul/sc_5/1/1.1/asses_object_1.1.rds"))
    unique_annotations <- unique(object@meta.data$main_new)
    #maxFractions <- numeric(length(unique_annotations)) 
for (i in 1:length(unique_annotations)) {
data_row <- as.matrix(table[i, ])
reshaped_data <- matrix(data_row, ncol = 10, byrow = FALSE)
result_matrix_reversed <- reshaped_data[, ncol(reshaped_data):1]
jpeg(paste0("/gpfs0/estiyl/users/tomul/sc_5/results/boxplot_sc_5_cluster_",i,".jpg"), width = 800, height = 600)
boxplot(result_matrix_reversed,
        main = paste0("Boxplot for cluster ", unique_annotations[i]),
        ylab = "Max Fraction",
        outline = FALSE,
        names = c("2879","6379","9879","13379", "16879", "20379", "23879","27379","30879","34379"),
        boxwex = 0.3,
        ylim = c(0, max(result_matrix_reversed, na.rm = TRUE)))
dev.off()


# making violin plot for each cluster
  data <- melt(result_matrix_reversed)
  colnames(data) <- c("Row", "Column", "Value")
  data$Column <- factor(data$Column, levels=1:10, labels=c("2879","6379","9879","13379", "16879", "20379", "23879","27379","30879","34379"))
  
jpeg(paste0("/gpfs0/estiyl/users/tomul/sc_5/results/violin_sc_5_cluster_", i, ".jpg"), width = 800, height = 600)

# Create the violin plot
p <- ggplot(data, aes(x=factor(Column), y=Value)) +
  geom_violin(trim=TRUE) +
  labs(x="Column", y="Value") +
  theme_minimal() +
  ggtitle(paste0("violin plot for cluster ", unique_annotations[i]))
  print(p)
  dev.off()


}





# Create an empty matrix to store the result
result_matrix <- matrix(NA, nrow = 90, ncol = 10)

for (i in seq(1, 90, by =  9)) {
  # Select the 10 columns
  selected_columns <- table[, i:(i + 8), drop = FALSE]
  
  # Bind the selected columns into the result matrix
  result_matrix[, (i %/%9 ) + 1] <- as.matrix(selected_columns)
}
result_matrix_reversed <- result_matrix[, ncol(result_matrix):1]
write.csv(result_matrix_reversed , paste0("/gpfs0/estiyl/users/tomul/sc_5/results/result_matrix_reversed.csv"))
#medians <- sapply(result_matrix_reversed, median, na.rm = TRUE)
#medians <- medians[2:length(medians)]
#cell_numbers <- c(460, 1460, 2460, 3460, 4460, 5460, 6460, 7460, 8460, 9460, 10460, 11460, 12460, 13460, 14460)
  # Perform linear regression using the first i data points
#lm_model <- lm(medians ~ cell_numbers)
#r_squared <- summary(lm_model)$r.squared
#r_squared_text <- paste0(" (R^2 = ", round(r_squared, digits=3), ")")
jpeg(paste0("/gpfs0/estiyl/users/tomul/sc_5/results/boxplot_sc_5_All_clusters.jpg"), width = 800, height = 600)
boxplot(result_matrix_reversed,
        main = paste0("Boxplot for all clusters"),
        ylab = "Max Fraction",
        outline = FALSE,
        names = c("2879","6379","9879","13379", "16879", "20379", "23879","27379","30879","34379"),
        las = 2,
        boxwex = 0.3,
        ylim = c(0, max(result_matrix_reversed, na.rm = TRUE)))
        

dev.off()

#x_values <- c(460, 1460, 2460, 3460, 4460, 5460, 6460, 7460, 8460, 9460, 10460, 11460, 12460, 13460, 14460)
#log_x_values <- log(x_values)
#jpeg(paste0("/gpfs0/estiyl/users/tomul/adipo10/results/boxplot_adipo10_All_clusters_Log.jpg"), width = 800, height = 600)
#boxplot(result_matrix_reversed,
#        main = paste0("Boxplot for all clusters"),
#        ylab = "Max Fraction",
#        outline = FALSE,
#        names = log_x_values,
#        las = 2,
#        boxwex = 0.3,
#        ylim = c(0, max(result_matrix_reversed, na.rm = TRUE)))
#dev.off()


data <- melt(result_matrix_reversed)
colnames(data) <- c("Row", "Column", "Value")
data$Column <- factor(data$Column, levels=1:10, labels=c("2879","6379","9879","13379", "16879", "20379", "23879","27379","30879","34379"))
  
jpeg(paste0("/gpfs0/estiyl/users/tomul/sc_5/results/violin_sc_5_all_clusters.jpg"), width = 800, height = 600)
# Create the violin plot
p <- ggplot(data, aes(x=factor(Column), y=Value)) +
  geom_violin(trim=TRUE) +
  labs(x="Column", y="Value") +
  theme_minimal() +
  ggtitle("violin plot for all clusters")
  print(p)
  dev.off()
  
  
  
#  ############################################################################################
#  selected_rows <- c(9, 7, 3, 16, 15, 12, 11)
#final_table_subset <- table[selected_rows,, drop=FALSE ]
#write.csv(final_table_subset , paste0("/storage16/projects/Noa/result/combinfinal_table_subset.csv"))
#  
#
## Create an empty matrix to store the result
#result_matrix <- matrix(NA, nrow = 80, ncol = 7)
#
## Loop through each set of 10 columns
#for (i in seq(1, 35, by =  5)) {
#  # Select the 10 columns
#  selected_columns <- final_table_subset[, i:(i + 4), drop = FALSE]
#  
#  # Bind the selected columns into the result matrix
#  result_matrix[, (i %/%5 ) + 1] <- as.matrix(selected_columns)
#}
#result_matrix_reversed <- result_matrix[, ncol(result_matrix):1]
#write.csv(result_matrix_reversed , paste0("/storage16/projects/Noa/result/result_matrix_reversed.csv"))
##medians <- sapply(result_matrix_reversed, median, na.rm = TRUE)
##medians <- medians[2:length(medians)]
##cell_numbers <- c(460, 1460, 2460, 3460, 4460, 5460, 6460, 7460, 8460, 9460, 10460, 11460, 12460, 13460, 14460)
#  # Perform linear regression using the first i data points
##lm_model <- lm(medians ~ cell_numbers)
##r_squared <- summary(lm_model)$r.squared
##r_squared_text <- paste0(" (R^2 = ", round(r_squared, digits=3), ")")
#jpeg(paste0("/storage16/projects/Noa/result/boxplot_sc_5_All_Adipo.jpg"), width = 800, height = 600)
#boxplot(result_matrix_reversed,
#        main = paste0("Boxplot for all Adipo"),
#        ylab = "Max Fraction",
#        outline = FALSE,
#        names = c("2000", "7000", "12000", "17000","22000","27000","32000"),
#        las = 2,
#        boxwex = 0.3,
#        ylim = c(0, max(result_matrix_reversed, na.rm = TRUE)))
#        
#
#dev.off()
#
##x_values <- c(460, 1460, 2460, 3460, 4460, 5460, 6460, 7460, 8460, 9460, 10460, 11460, 12460, 13460, 14460)
##log_x_values <- log(x_values)
##jpeg(paste0("/gpfs0/estiyl/users/tomul/adipo10/results/boxplot_adipo10_All_clusters_Log.jpg"), width = 800, height = 600)
##boxplot(result_matrix_reversed,
##        main = paste0("Boxplot for all clusters"),
##        ylab = "Max Fraction",
##        outline = FALSE,
##        names = log_x_values,
##        las = 2,
##        boxwex = 0.3,
##        ylim = c(0, max(result_matrix_reversed, na.rm = TRUE)))
##dev.off()
#
#
#data <- melt(result_matrix_reversed)
#colnames(data) <- c("Row", "Column", "Value")
#data$Column <- factor(data$Column, levels=1:7, labels=c( "2000", "7000", "12000", "17000","22000","27000","32000"))
#  
#jpeg(paste0("/storage16/projects/Noa/result/violin_sc_5_all_Adipo.jpg"), width = 800, height = 600)
## Create the violin plot
#p <- ggplot(data, aes(x=factor(Column), y=Value)) +
#  geom_violin(trim=TRUE) +
#  labs(x="Column", y="Value") +
#  theme_minimal() +
#  ggtitle("violin plot for all Adipo")
#  print(p)
#  dev.off()