# Original Source: Adamic & Glance (2005), The political blogosphere and the 2004 U.S. election: divided they blog, Proceedings of the 3rd international workshop on Link discovery.
# txt File Source: https://github.com/p-wolfe/network-histogram-code/blob/master/polblogs.txt

polblog_edgelist <- read.table("data-raw/polblogs.txt")
polblog <- igraph::graph_from_edgelist(as.matrix(polblog_edgelist), directed = FALSE)
polblog <- igraph::delete_vertices(polblog, igraph::degree(polblog) == 0) #remove isolated nodes

usethis::use_data(polblog, polblog_edgelist, overwrite = TRUE)