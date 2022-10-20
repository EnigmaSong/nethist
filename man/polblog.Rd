\name{polblog}
\alias{polblog}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Political blog 
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Political blog data set from Adamic and Glance (2005).
}
\usage{data("polblog")}
\format{
  A list of edges.
  \describe{
    \item{\code{V1}}{a numeric vector of tail vertices.}
    \item{\code{V2}}{a numeric vector of head vertices.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The adjancey matrix is sparse and assymetric, so it needs to be symmetrized if a method requires a symmetric matrix.
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
data(polblog)
## maybe str(polblog) ; plot(polblog) ...
}
\keyword{datasets}
