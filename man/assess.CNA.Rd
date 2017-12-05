\name{assess.CNA}
\alias{assess.CNA}
\title{
Significance assessment of copy number alterations
}
\description{
Starting from a matrix of coverages, the corresponding matrix of copy numbers is calculated.
A null model for significance assessment of copy number alterations is fitted. Each amplicon
in each sample is assessed for significance. Summarized copy numbers and p-values for genes
are calculated as descibed below.
}
\usage{
assess.CNA(coverage.target, coverage.source=NULL, method.pooled="amplicon", thres.cov=100)
}
\arguments{
  \item{coverage.target}{
A numeric matrix containing the target coverages of each amplicon (rows) in each sample (columns).
The target data are investigated for copy number alterations.
}
  \item{coverage.source}{
A numeric matrix containing the source coverages of each amplicon (rows) in each sample (columns).
The source data are used to fit a null model. If \code{NULL}, the target data are used to fit the
null model.
}
  \item{method.pooled}{
Method used for the estimation of the null model. Either one common null model for all amplicons
(\code{pooled}) or individual null models for each of the amplicons (\code{amplicon}) are fitted.
}
  \item{thres.cov}{
Theshold for the minimal mean coverage of an amplicon to be included in the analysis.
}

}
\value{
List containing the following elements: Matrix of copy numbers with the estimated null model (\code{"model"}), estimates of copy numbers (\code{"CN.a"} and \code{"CN.g"}) for amplicons and genes as well as p-values of copy number alterations (\code{"P.a"} and \code{"P.g"}) for amplicons and genes. Copy numbers for genes are calculated as average of the copy numbers of all amplicons interrogating the gene, p-values for genes are calculated using Fisher's method.
}
\examples{
\dontrun{
data(coverage)
CNA <- assess.CNA(coverage)
}
\dontshow{
data(coverage)
CN <- assess.CNA(coverage[1:3, ])
}
}
