\name{call.amplicons}
\alias{call.amplicons}
\title{
Amplicon-wise call of copy number alterations
}
\description{
Call of significant copy number alterations including different methods to correct
for multiple hypothesis testing
}
\usage{
call.amplicons(CN, direction = "gain", method.p = "p_samples", thres.p = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CN}{
Matrix of copy numbers with attributes as decribed in \code{\link{calculate.CN}}.
}
  \item{direction}{
Detection of gene amplifications ("gain") or gene deletions ("loss").
}
  \item{method.p}{
Usage of uncorrected p-values ("p"), p-values corrected for samples ("p_samples")
or p-values corrected for samples and amplicons ("p_samples_amplicons").
}
  \item{thres.p}{
Significance level for calling of copy number alterations.
}
}
\value{
Indicator matrix with the detected copy number alterations in each amplicon ("amplicon")
and in each sample ("sample") as attributes.
}
\examples{
\dontrun{
data(coverage)
CN <- calculate.CN(coverage)
calls <- call.amplicons(CN, direction="gain")
}
\dontshow{
data(coverage)
CN <- calculate.CN(coverage[1:3, ])
calls <- call.amplicons(CN, direction="gain")
}
}
