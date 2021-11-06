# Find DEGs between old and young cell clusters
getDEGs <- function(seurat.object = seurat.object, ident.1 = ident.1, ident.2 = ident.2){
  
  n <- levels(Idents(seurat.object))
  Idents(seurat.object) <- "group"
  DefaultAssay(seurat.object) <- "RNA"
  
  DEGs <- list()
  for (i in n) {
    print(paste("Finding DEGs of cluster ",i))
    DEGs[[i]] <- FindMarkers(seurat.object, ident.1 = paste0(i,"_",ident.1), ident.2 = paste0(i,"_",ident.2), logfc.threshold = 0.25, verbose = FALSE)
    Gene <- as.character(rownames(DEGs[[i]]))
    DEGs[[i]] <- cbind(Gene,DEGs[[i]])
    DEGs[[i]] <- DEGs[[i]][order(DEGs[[i]]$avg_logFC,decreasing = T),]
    DEGs[[i]] <- DEGs[[i]][DEGs[[i]]$p_val_adj < 0.05, ]
  }
  
  return(DEGs)
  
}



# Get common up-regulated genes in old cell clusters
getCommonUpDEGs <- function(DEG_list){
  
  list_of_DEGs <- list()
  for (i in 1:length(DEG_list)) {
    DEG_data <- as.data.frame(DEG_list[[i]])
    list_of_DEGs[[i]] <- rownames(DEG_data)[which(DEG_data$avg_logFC>0)]
  }
  
  num_dec <- list()
  for (i in 1:length(list_of_DEGs)) {
    num_dec[[i]] <- as.numeric(rep(0,length(list_of_DEGs[[i]])))
  }
  
  for (i in 1:length(list_of_DEGs)) {
    DEGs <- list_of_DEGs[[i]]
    for (j in 1:length(DEGs)) {
      for (k in 1:length(list_of_DEGs)) {
        if (DEGs[[j]] %in% list_of_DEGs[[k]]) {
          num_dec[[i]][j] <- num_dec[[i]][j]+1
        }
      }
    }
  }
  
  common_DEGs <- list()
  for (i in 1:length(list_of_DEGs)) {
    common_DEGs[[i]] <- list_of_DEGs[[i]][which(num_dec[[i]]>=4)]
  }
  common_DEGs <- Reduce(union, common_DEGs)
  return(common_DEGs)
}



# Get common down-regulated genes in old cell clusters
getCommonDownDEGs <- function(DEG_list){
  
  list_of_DEGs <- list()
  for (i in 1:length(DEG_list)) {
    DEG_data <- as.data.frame(DEG_list[[i]])
    list_of_DEGs[[i]] <- rownames(DEG_data)[which(DEG_data$avg_logFC<0)]
  }
  
  num_dec <- list()
  for (i in 1:length(list_of_DEGs)) {
    num_dec[[i]] <- as.numeric(rep(0,length(list_of_DEGs[[i]])))
  }
  
  for (i in 1:length(list_of_DEGs)) {
    DEGs <- list_of_DEGs[[i]]
    for (j in 1:length(DEGs)) {
      for (k in 1:length(list_of_DEGs)) {
        if (DEGs[[j]] %in% list_of_DEGs[[k]]) {
          num_dec[[i]][j] <- num_dec[[i]][j]+1
        }
      }
    }
  }
  
  common_DEGs <- list()
  for (i in 1:length(list_of_DEGs)) {
    common_DEGs[[i]] <- list_of_DEGs[[i]][which(num_dec[[i]]>=4)]
  }
  common_DEGs <- Reduce(union, common_DEGs)
  return(common_DEGs)
}



# Get unique up-regulated genes in old cell clusters
getUniqueUpGenes <- function(DEG_list){
  
  up_DEGs <- list()
  for (i in 1:length(DEG_list)) {
    up_DEGs[[i]] <- rownames(as.data.frame(DEG_list[[i]]))[which(as.data.frame(DEG_list[[i]])$avg_logFC>0)]
  }
  
  
  list_of_DEGs <- list()
  for (i in 1:length(up_DEGs)){
    list_of_DEGs[[i]] <- setdiff(up_DEGs[[i]], Reduce(union, up_DEGs[-i]))
  }
  
  return(list_of_DEGs)
}



# Get unique down-regulated genes in old cell clusters
getUniqueDownGenes <- function(DEG_list){
  
  down_DEGs <- list()
  for (i in 1:length(DEG_list)) {
    down_DEGs[[i]] <- rownames(as.data.frame(DEG_list[[i]]))[which(as.data.frame(DEG_list[[i]])$avg_logFC<0)]
  }
  
  
  list_of_DEGs <- list()
  for (i in 1:length(down_DEGs)){
    list_of_DEGs[[i]] <- setdiff(down_DEGs[[i]], Reduce(union, down_DEGs[-i]))
  }
  
  return(list_of_DEGs)
}




