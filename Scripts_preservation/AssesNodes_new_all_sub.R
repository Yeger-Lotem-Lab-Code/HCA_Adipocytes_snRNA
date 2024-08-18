suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(ranger))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

AssessSplit <- function(
  object,
  node,
  cluster1,
  cluster2,
  genes.training = NULL,
  print.output = TRUE,
  ...
) {
  genes.training <- SetIfNull(x = genes.training, default = rownames(x = object@assays$RNA@data))
  genes.training <- intersect(x = genes.training, rownames(x = object@assays$RNA@data))
  if (!length(x = genes.training)) {
    stop("None of the genes provided are in the data")
  }
  tree <- object@tools$BuildClusterTree
  if (!missing(x = node)) {
    if (!missing(x = cluster1) || !missing(x = cluster2)) {
      warning("Both node and cluster IDs provided. Defaulting to using node ID")
    }
    possible.nodes <- c(
      DFT(tree = tree, node = tree$edge[,1][1]),
      tree$edge[,1][1]
    )
    if (!node %in% possible.nodes) {
      stop("Not a valid node")
    }
    split <- tree$edge[which(x = tree$edge[,1] == node), ][,2]
    group1 <- DFT(tree = tree, node = split[1], only.children = TRUE)
    group2 <- DFT(tree = tree, node = split[2], only.children = TRUE)
    if (any(is.na(x = group1))) {
      group1 <- split[1]
    }
    if (any(is.na(x = group2))) {
      group2 <- split[2]
    }
  } else {
    group1 <- cluster1
    group2 <- cluster2
  }
  group1.cells <- WhichCells(object = object, ident = group1)
  group2.cells <- WhichCells(object = object, ident = group2)
  assess.data <- subset(
    x = object,
    cells = c(group1.cells, group2.cells)
  )
  assess.data <- SetIdent(
    object = assess.data,
    cells = group1.cells,
    value = "g1"
  )
  assess.data <- SetIdent(
    object = assess.data,
    cells = group2.cells,
    value = "g2"
  )
  rfc <- BuildRFClassifier(
    object = assess.data,
    # training.genes = assess.data@var.genes,
    training.genes = genes.training,
    training.classes = assess.data@active.ident,
    ...
  )
  oobe <- rfc$prediction.error
  if (print.output) {
    message(paste0("Out of Bag Error: ", round(x = oobe, digits = 4) * 100, "%"))
  }
  return(oobe)
}
AssessNodes <- function(
  object,
  node.list,
  all.below = FALSE,
  genes.training = NULL
) {
  genes.training <- SetIfNull(x = genes.training, default = rownames(x = object@assays$RNA@data))
  genes.training <- intersect(x = genes.training, rownames(x = object@assays$RNA@data))
  if (!length(x = genes.training)) {
    stop("None of the genes provided are in the data")
  }
  tree <- object@tools$BuildClusterTree
  if (missing(x = node.list)) {
    node.list <- GetAllInternalNodes(tree = tree)
  } else {
    possible.nodes <- GetAllInternalNodes(tree = tree)
    if (any(!node.list %in% possible.nodes)) {
      stop(paste(
        node.list[!(node.list %in% possible.nodes)],
        "not valid internal nodes"
      ))
    }
    if (length(x = node.list == 1) && all.below) {
      node.list <- c(node.list, DFT(tree = tree, node = node.list))
    }
  }
  oobe <- pbsapply(
    X = node.list,
    FUN = function(x) {
      return(AssessSplit(
        object = object,
        node = x,
        genes.training = genes.training,
        print.output = FALSE,
        verbose = FALSE
      ))
    }
  )
  return(data.frame(node = node.list, oobe))
}

SetIfNull <- function(x, default) {
  if (is.null(x = x)) {
    return(default)
  } else {
    return(x)
  }
}

GetAllInternalNodes <- function(tree) {
  return(c(tree$edge[1, 1], DFT(tree = tree, node = tree$edge[1, 1])))
}

DFT <- function(
  tree,
  node,
  path = NULL,
  include.children = FALSE,
  only.children = FALSE
) {
  if (only.children) {
    include.children = TRUE
  }
  children <- which(x = tree$edge[, 1] == node)
  child1 <- tree$edge[children[1], 2]
  child2 <- tree$edge[children[2], 2]
  if (child1 %in% tree$edge[, 1]) {
    if(! only.children){
      path <- c(path, child1)
    }
    path <- DFT(
      tree = tree,
      node = child1,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <-c(path, child1)
    }
  }
  if (child2 %in% tree$edge[, 1]) {
    if (! only.children) {
      path <- c(path, child2)
    }
    path <- DFT(
      tree = tree,
      node = child2,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <- c(path, child2)
    }
  }
  return(path)
}

BuildRFClassifier <- function(
  object,
  training.genes = NULL,
  training.classes = NULL,
  verbose = TRUE,
  ...
) {
  PackageCheck('ranger')
  training.classes <- as.vector(x = training.classes)
  training.genes <- SetIfNull(
    x = training.genes,
    default = rownames(x = object@assays$RNA@data)
  )
  training.data <- as.data.frame(
    x = as.matrix(
      x = t(
        x = as.matrix(object@assays$RNA@data[training.genes, ])
      )
    )
  )
  training.data$class <- factor(x = training.classes)
  if (verbose) {
    message("Training Classifier ...")
  }
  classifier <- ranger::ranger(
    data = training.data,
    dependent.variable.name = "class",
    classification = TRUE,
    write.forest = TRUE,
    ...
  )
  return(classifier)
}

PackageCheck <- function(..., error = TRUE) {
  pkgs <- unlist(x = c(...), use.names = FALSE)
  package.installed <- vapply(
    X = pkgs,
    FUN = requireNamespace,
    FUN.VALUE = logical(length = 1L),
    quietly = TRUE
  )
  if (error && any(!package.installed)) {
    stop(
      "Cannot find the following packages: ",
      paste(pkgs[!package.installed], collapse = ', '),
      ". Please install"
    )
  }
  invisible(x = package.installed)
}

MergeNode <- function(object, node.use, rebuild.tree = FALSE, ...) {
  object.tree <- object@tools$BuildClusterTree
  node.children <- DFT(
    tree = object.tree,
    node = node.use,
    include.children = TRUE
  )
  node.children <- intersect(x = node.children, y = levels(x = object@active.ident))
  children.cells <- WhichCells(object = object, ident = node.children)
  if (length(x = children.cells > 0)) {
    object <- SetIdent(
      object = object,
      cells = children.cells,
      value = min(node.children)
    )
  }
  if (rebuild.tree) {
    object <- BuildClusterTree(object = object, ...)
  }
  return(object)
}
 
# part one
# input_path <- "G:/My Drive/MSc-Phd/adipo/adipo_fraction_analysis/compare_adipocytes/"
# output_path <- paste0(input_path, "assessed/")
for (i in 6:6) {
  for (j in 1:1) {
  
output_path <- paste0("/gpfs0/estiyl/users/tomul/sc_5/",i,"/",i,".",j,"/")
file_path <- paste0("/gpfs0/estiyl/users/tomul/sc_5/",i, "/",i, ".",j, "/harmony_",i, ".",j,".rds")
asses_object <- readRDS("/gpfs0/estiyl/users/tomul/sc_5/6/6.1/harmony_6.1.rds")
asses_object <- BuildClusterTree(asses_object,reorder = T, reorder.numeric = T, dims = 1:50) 
PlotClusterTree(asses_object)
# need to change to umap_orig?
DimPlot(asses_object, reduction = "umap_orig", label = T, pt.size = 1.5, label.size = 8)+theme(legend.position = "none")
node.scores <- AssessNodes(asses_object, genes.training = VariableFeatures(asses_object))
node.scores <- node.scores[order(node.scores$oobe,decreasing = T),]
write.csv(node.scores, file = paste(output_path, "asses_object_node_",i,".",j,".scores.csv", sep=""))
  
# part two

# dir.create(output_path)
print(node.scores)
num_rows_greater_than_0.05 <- nrow(node.scores[node.scores$oobe > 0.05, ])
nodes.merge <- node.scores[1:num_rows_greater_than_0.05,] #choose insignificant splits (x)
print(nodes.merge)
nodes.to.merge <- sort(nodes.merge$node, decreasing = T)
print(nodes.to.merge)
for (n in nodes.to.merge){
  asses_object <- MergeNode(asses_object, n)}
  idents <- levels(x=asses_object)
if(length(idents)>1){
asses_object <- BuildClusterTree(asses_object,reorder = T, reorder.numeric = T, dims = 1:50) 
PlotClusterTree(asses_object)

node.scores <- AssessNodes(asses_object, genes.training = VariableFeatures(asses_object))
node.scores <- node.scores[order(node.scores$oobe,decreasing = T),]
write.csv(node.scores, file = paste(output_path, "asses_object_node_",i,".",j,".scores2.csv", sep=""))
rm(node.scores, nodes.merge,nodes.to.merge)
}
p2 <- DimPlot(asses_object, reduction = "umap_orig", label = T, pt.size = 1, label.size = 8) +
  theme(legend.position = "none")
jpeg(paste(output_path, "asses_object_",i,".",j,"_umap_byCluster2.jpeg", sep = ""), width = 1500, height = 1500)
print(p2)
dev.off()
saveRDS(asses_object, file = paste(output_path, "asses_object_",i,".",j,".rds", sep = ""))

  }
}

