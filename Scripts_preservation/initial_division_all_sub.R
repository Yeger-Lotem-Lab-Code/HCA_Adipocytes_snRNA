object <- readRDS("/gpfs0/estiyl/users/lazaresc/figs4paper/sc5_diet.rds")
object@reductions$umap_orig <- object@reductions$umap

# barcode names
col_names <- colnames(object@assays$RNA@counts)

# random shuffle of barcodes vector
shuffled_names <- sample(col_names)
vec_list <- list()
size_vec <- length(col_names) # 37879 in sc5_diet
sub_vector <- shuffled_names[3501:size_vec, drop = FALSE] # taking all cells besides the first 5000
vec_list[[1]] <- sub_vector
sub_object <- subset(object,cells = sub_vector)
saveRDS(sub_object, file = paste0("/storage16/projects/Tom_Noa/sc_5/1/1.1/sc_1.1.rds") )

# saving all the objects with 3500 less cells that were randomly chosen.
for (i in 2:10) {
  start_index <- (i-1)*3500
  end_index <- (i*3500+1)
  selected_cols <- c(1:start_index , end_index:size_vec)
  sub_vector_I <- shuffled_names[selected_cols, drop = FALSE]
  vec_list[[i]] <- sub_vector_I
  sub_object <- subset(object,cells = sub_vector_I)
  file_name <- paste0("/storage16/projects/Tom_Noa/sc_5/", i, "/", i, ".1/sc_", i, ".1.rds")
  saveRDS(sub_object, file = file_name)
}


for (i in 1:10) {
  sub_vector_I <- vec_list[[i]]
  size_vec <- length(sub_vector_I)
  subsetSize = size_vec
  for (j in 2:10){
    

    subsetSize = subsetSize - 3500
    
    sub_vector_I <- sample(sub_vector_I, size = subsetSize, replace=F)
    sub_object <- subset(object,cells = sub_vector_I)
    file_name <- paste0("/storage16/projects/Tom_Noa/sc_5/", i, "/", i, ".", j, "/sc_", i, ".", j, ".rds")
    saveRDS(sub_object, file = file_name)
      
  }
}