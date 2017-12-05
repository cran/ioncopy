\name{call.CNA}
\alias{call.CNA}
\title{
Gene-wise or amplicon-wise calling of copy number alterations
}
\description{
Copy number alteration calling after (possible) multiple testing correction for samples and/or amplicons/genes.
P-values are controlled either for the number of samples, the number of amplicons/genes or both.
Different methods of multiple testing can be chosen to control either FWER or FDR.
}
\usage{
call.CNA(CNA, analysis.mode="gene-wise", method.p="samples_genes/amplicons",
method.mt="Bonferroni", thres.p=0.05, sig.call=0)
}
\arguments{
  \item{CNA}{
List of CNA assessments generated by \code{\link{assess.CNA}}.
}
  \item{analysis.mode}{
The Mode of the analysis: (\code{"gene-wise"}) or (\code{"amplicon-wise"})
}
  \item{method.p}{
The multiple testing method used for detection:
Usage of uncorrected p-values (\code{"p"}), p-values corrected for samples (\code{"p_samples"}),
p-values corrected for amplicons/genes (\code{"p_genes/amplicons"}) or p-values corrected for samples and
amplicons/genes (\code{"p_samples_genes/amplicons"}).
}
  \item{method.mt}{
Method for multiple testing correction: must be equal to \code{Bonferroni}, \code{Holm},
\code{Hochberg}, \code{SidakSS}, \code{SidakSD}, \code{BH}, \code{BY}, \code{ABH} or \code{TSBH}
as described in the package \code{multtest}.
}
  \item{thres.p}{
Significance level for calling of copy number alterations.
}
  \item{sig.call}{
An integer >= 0. Only used if \code{analysis.mode="gene-wise"}. Required minimum number of amplicons supporting a gene CNA call.
}

}

\value{
Table containing the status (GAIN, LOSS or NORMAL) of each amplicon or gene in each sample (\code{"tab"}). Matrix of copy number estimates in each sample and each amplicon or gene (\code{"CN"}). Indicator matrix of detected gains in each sample and each amplicon or gene (\code{"gain"}) and the same for losses (\code{"loss"}).
}
\examples{
\dontrun{
data(coverage)
CN <- calculate.CN(coverage)
CNA <- assess(CN)
calls <- call.CNA(CNA)
}
\dontshow{
data(coverage)
CN <- assess.CNA(coverage)
CNA <- call.CNA(CN, analysis.mode = "amplicon-wise", method.p = "samples", method.mt = "BH", thres.p = 0.05, sig.call= 1)
}
}