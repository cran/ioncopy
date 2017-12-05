\name{calculate.CN}
\alias{calculate.CN}
\title{
Calculation of copy numbers
}
\description{
Starting from a matrix of coverages, the corresponding matrix of copy numbers is calculated.
The caluculation includes a sample normalization step and a amplicon normalization step.
}
\usage{
calculate.CN(coverage, scale.amplicon=NULL)
}
\arguments{
  \item{coverage}{
A numeric matrix containing the coverages of each amplicon (rows) in each sample (columns).
}
  \item{scale.amplicon}{
Method for amplicon normalization. If \code{NULL} the scaling factors are estimated from \code{coverage}.
}
}
\value{
Matrix of copy numbers with the amplicon scaling factors (\code{"scale.amplicon"}) as attribute.
}
\examples{
\dontrun{
data(coverage)
CN <- calculate.CN(coverage)
}
\dontshow{
data(coverage)
CN <- calculate.CN(coverage[1:3, ])
}
}
