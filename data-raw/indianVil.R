# Original Source: https://doi.org/10.7910/DVN/U3BIHX
# From the source, we only use Village No. 40.
files_adj_mats <- list.files(path = "data-raw", full.names = TRUE)
layer_names <- c("Borrow money",
                 "give advice",
                 "help decision",
                 "kero rice come",
                 "kero rice go",
                 "lend money",
                 "medic",
                 "nonrel",
                 "rel",
                 "temple company",
                 "visit come",
                 "visit go")
files_HH_networks <- unlist(sapply(files_adj_mats, grep, 
                            pattern = paste0("_HH_vilno_40.csv"), value = TRUE),
                            use.name = FALSE)
#Load csv files (adj_borrowmoney_HH_vilno_40.csv, ....)
nets <- vector("list", length(files_HH_networks))
for(i in 1:length(files_HH_networks)){
  print(files_HH_networks[i])
  loaded_adj_mat <- as.matrix(read.csv(files_HH_networks[i], header = FALSE))
  nets[[i]] <- igraph::graph_from_adjacency_matrix(loaded_adj_mat,
                                                mode = "max", diag = FALSE)
}
#Covert list to array
n_nodes <- igraph::vcount(nets[[1]])
IndianVil <- array(NA, dim = c(n_nodes,n_nodes,12))
for(l in 1:12){
  IndianVil[,,l] <- igraph::as_adjacency_matrix(nets[[l]], sparse=FALSE)
}
#Remove nodes with zero degree across all layers.
glob_deg_net_array <- rowSums(IndianVil)
zero_deg_nodes <- which(glob_deg_net_array==0)
IndianVil <- IndianVil[-zero_deg_nodes,-zero_deg_nodes,]
dimnames(IndianVil)[[3]] <- layer_names

usethis::use_data(IndianVil, overwrite = TRUE)

