object <-  readRDS(file = paste0("/gpfs0/estiyl/users/tomul/sc_5/1/1.1/asses_object_1.1.rds"))
unique_annotations <- unique(object@meta.data$main_new)

for (n in 1:10) {
  for (m in 1:10) {
    maxFractions <- numeric(length(unique_annotations))
    table <- read.csv(file = paste0("/gpfs0/estiyl/users/tomul/sc_5/", n, "/", n, ".", m, "/output.csv"))
    
    # looping thru all clusters
    for (i in 1:length(unique_annotations)) {
    
    
      max_value <- max(table$table2, na.rm = TRUE)# max_value represents number of clusters we have in subgroup i.j . if we found 6 clusters then max_value = 6
      maxFraction <- 0
      origClusterName <- unique_annotations[i]
      for (j in 1:max_value) {
        numOfCells <- sum(table$table1 == origClusterName & table$table2 == j, na.rm = TRUE)
        totalNum <- sum(table$table2 == j, na.rm = TRUE)
        maxFraction <- max(maxFraction, numOfCells / totalNum)
      }
      maxFractions[i] <- maxFraction
    }
    
    maxFractions <- data.frame(maxFractions)
    write.csv(maxFractions, paste0("/gpfs0/estiyl/users/tomul/sc_5/", n, "/", n, ".", m, "/maxFractions.csv"))
  }
}