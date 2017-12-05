\name{read.coverages}
\alias{read.coverages}
\title{
Helper function to read coverages from file system
}
\description{
Coverages are read from tab separated files and stored in a coverage matrix.
All coverage files need to refer to the same panel of targets.
}
\usage{
read.coverages(chip.names, file.names, anno.col="Target")
}
\arguments{
  \item{chip.names}{
A character vector comprizing the chip names.
}
  \item{file.names}{
A character vector compizing the names of tab separated files.
Each file should contain the coverages of the same sequencing
panel (rows) in a cohort of samples (columns).
}
  \item{anno.col}{
The column in input files that compizes the names of the targets.
The targets must belong to the same sequencing panel for all input files.
}
} 
\value{
Matrix of coverages of each target sequence (rows) in each sample (colums).
}
