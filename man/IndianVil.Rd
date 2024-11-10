\name{IndianVil}
\alias{IndianVil}
\docType{data}
\title{
  %%   ~~ data name/kind ... ~~
    Indian Village socioeconomic networks
}
\description{
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
 A processed version of the Indian village dataset from Banerjee et al. (2013), representing socio-economic networks at the individual level for village ID 40.
}
\usage{
  data(IndianVil)
}
\format{
  \code{IndianVil} is an array with size of 231 x 231 x 12.
}
\details{
  %%  ~~ If necessary, more details than the __description__ above ~~
  The array represents the socio-economic relationships of individuals from Village ID 40. Nodes with zero degrees across all layers were removed, reducing the total number of nodes from 241 to 231.
}
\source{
  %%  ~~ reference to a publication or URL from which the data were obtained ~~
  https://doi.org/10.7910/DVN/U3BIHX
}
\references{
  %%  ~~ possibly secondary sources and usages ~~
    Banerjee, A., Chandrasekhar, A. G., Duflo, E., & Jackson, M. O. (2013). The diffusion of microfinance. Science, 341(6144), 1236498.
    Song, Y. & Olhede, S.C. (2024+) 
}
\examples{
  data(InvianVil)
  
  InvianVil
}
\keyword{datasets}
