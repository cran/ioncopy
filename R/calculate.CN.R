calculate.CN <-
function(coverage, scale.amplicon=NULL) {
  scale.sample <- apply(coverage, 2, median)
  C <- scale(coverage, center=FALSE, scale=scale.sample/mean(scale.sample))
  if (is.null(scale.amplicon)) scale.amplicon <- apply(C, 1, median) / 2
  CN <- t(scale(t(C), center=FALSE, scale=scale.amplicon))
  attr(CN, 'scale.amplicon') <- scale.amplicon
  attr(CN, 'scale.sample') <- scale.sample
  return(CN)
}
