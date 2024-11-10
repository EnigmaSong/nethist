\name{polblog}
\alias{polblog}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Political blog 
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
A pre-processed version of the Political Blog dataset from Olhede and Wolfe (2014), based on the original data by Adamic and Glance (2005). The provided igraph object represents the graph after removing zero-degree nodes and simplifying the graph structure. The graph is derived from the edge list \code{polblog_edgelist}.
}
\usage{
data(polblog)
}
\format{
\code{polblog} is an igraph object, which is pre-processed data symmetrize and simplified without self-loops.
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The adjancency matrix from \code{polblog_edgelist} is sparse, assymetric, so it needs to be symmetrized if a method requires a symmetric matrix. Nodes with zero degree need to be removed. All edges are considered as undirected edges. Then, we get the undirected version of dataset as \code{polblog}.

First 586 blogs are liberal; remaining 638 are conservative.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
http://www-personal.umich.edu/~mejn/netdata/

https://github.com/p-wolfe/network-histogram-code
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Adamic, L. A., & Glance, N. (2005). The political blogosphere and the 2004 US election: divided they blog. In Proceedings of the 3rd international workshop on Link discovery (pp. 36-43)
}
\examples{
data(polblog_edgelist)
data(polblog)

#From polblog_edgelist to polblog
G<- igraph::graph_from_edgelist(as.matrix(polblog_edgelist), directed = FALSE)
G<- igraph::delete.vertices(G, igraph::degree(G) == 0)
G<- igraph::simplify(G)

polblog
}
\keyword{datasets}
