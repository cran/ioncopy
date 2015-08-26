\name{calculate.CN}
\alias{calculate.CN}
\title{
Calculation of copy numbers
}
\description{
Starting from a matrix of coverages, the corresponding matrix of copy numbers is calculated.
A null model for significance assessment of copy number alterations is fitted. Each amplicon
in each sample is assessed for significance. P-values are reported in three ways:
without multiple testing correction, correted for the number of samples and corrected for
both, the number of samples and the number of amplicons.
}
\usage{
calculate.CN(coverage, method.pooled = "amplicon", method.mt = "Bonferroni",
thres.cov = 100, thres.p = 0.05)
}
\arguments{
  \item{coverage}{
A numeric matrix containing the coverages of each amlicon (rows) in each of samples (columns).
}
  \item{method.pooled}{
Method used for the estimation of the null model. Either one common null model for all amplicons
("pooled") or individual null models for each of the amplicons ("amplicon") are fitted.
("amplicon")
}
  \item{method.mt}{
Method for multiple testing correction: must be equal to "Bonferroni", "Holm", "Hochberg",
"SidakSS", "SidakSD", "BH", "BY", "ABH" or "TSBH" as described in the package multtest.
}
  \item{thres.cov}{
Theshold for the minimal mean coverage of an amplicon to be included in the analysis.
}
  \item{thres.p}{
Significance level for calling copy number alterations.
}
}
\value{
Matrix of copy numbers with the estimated null model ("model") and the p-values
for each amplicon in each sample ("P") as attributes.
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
