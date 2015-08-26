\name{call.genes}
\alias{call.genes}
\title{
Gene-wise calling of copy number alterations
}
\description{
Gene calls are calculated as summary of the calls of all amplicons in the gene.
}
\usage{
call.genes(CN, direction = "gain", method.p.det = "p_samples_amplicons",
method.p.val = "p_samples", n.validated = 1, thres.p = 0.05)
}
\arguments{
  \item{CN}{
Matrix of copy numbers with attributes as decribed in \code{\link{calculate.CN}}.
}
  \item{direction}{
Detection of gene amplifications ("gain") or gene deletions ("loss").
}
  \item{method.p.det}{
Multiple testing method for detection: Usage of uncorrected p-values ("p"),
p-values corrected for samples ("p_samples")
or p-values corrected for samples and amplicons ("p_samples_amplicons").
}
  \item{method.p.val}{
Multiple testing method for validation: Usage of uncorrected p-values ("p"),
p-values corrected for samples ("p_samples")
or p-values corrected for samples and amplicons ("p_samples_amplicons").
}
  \item{n.validated}{
Number of required calls by additional amplicons in the same gene to consider
a detected copy number alteration as validated. 
}
  \item{thres.p}{
Significance level for calling of opy number alterations.
}
}
\value{
Indicator matrix with the detected copy number alterations in each gene ("gene")
and in each sample ("sample") as attributes.
}
\examples{
\dontrun{
data(coverage)
CN <- calculate.CN(coverage)
calls <- call.genes(CN, direction="gain")
}
\dontshow{
data(coverage)
CN <- calculate.CN(coverage[1:3, ])
calls <- call.genes(CN, direction="gain")
}
}
