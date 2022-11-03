\name{polblog}
\alias{polblog}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Political blog 
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Political blog data set from Olhede and Wolfe (2014), which is part of the original data set from Adamic and Glance (2005). The provided edgelist is from GitHub repository of Prof. Patrick Wolfe. The provided igraph object is the zero-degree removed and simplified graph from the edgelist. 

}
\usage{
data(polblog)
data(polblog_edgelist)
}
\format{
- polblog is an igraph object of the political blog dataset, which is pre-processed data symmetrize and simplified without self-loops.

-  polblog_edgelist is a matrix for the edge list from https://github.com/p-wolfe/network-histogram-code.
  \describe{
    \item{\code{V1}}{a numeric vector of tail vertices.}
    \item{\code{V2}}{a numeric vector of head vertices.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The adjancency matrix from polblog_edgelist is sparse, assymetric, so it needs to be symmetrized if a method requires a symmetric matrix.
First 586 blogs are liberal; remaining 638 are conservative. Nodes with zero degree need to be removed.
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

G<- igraph::graph_from_edgelist(as.matrix(polblog_edgelist), directed = FALSE)
G<- delete.vertice(G, degree(G) == 0)
G<- igraph::simplify(G)

polblog
}
\keyword{datasets}
