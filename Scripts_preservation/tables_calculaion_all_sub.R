# this script gives a table for each cell that tells the original cluster and the new cluster.
for (i in 1:10) {
  for (j in 1:10) {
object <-  readRDS(file = paste0("/gpfs0/estiyl/users/tomul/vis_10/",i,"/",i,".",j,"/asses_object_", i, ".", j, ".rds"))
table1 <- object@meta.data$main_new # need to get the original cluster here
table1 <- data.frame(table1)
table2 <- object@active.ident
table2 <- data.frame(table2)
table <- cbind(table1,table2)
write.csv(table , paste0("/gpfs0/estiyl/users/tomul/vis_10/",i,"/",i,".",j,"/output.csv"))
  }
}