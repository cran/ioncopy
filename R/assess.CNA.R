
assess.CNA <-
function(coverage.target, coverage.source=NULL, method.pooled="amplicon", thres.cov=100) {

  if (is.null(coverage.source)) coverage.source <- coverage.target
  CN.source <- calculate.CN(coverage.source)
  CN.target <- calculate.CN(coverage.target, scale.amplicon=attr(CN.source, "scale.amplicon"))

  scale.amplicons <- attr(CN.target, 'scale.amplicon')
  scale.samples <- attr(CN.target, 'scale.sample')

  m.coverage.amplicon <- apply(coverage.source, 1, mean)
  amplicons <- names(m.coverage.amplicon)[which(m.coverage.amplicon >= thres.cov)]
  namplicons <- length(amplicons)
  samples <- colnames(CN.target)
  nsamples <- length(samples)
  CN.source <- CN.source[amplicons, ]
  CN.target <- CN.target[amplicons, ]

 #1: MODEL MATRIX
  model <- matrix(nrow=namplicons+1, ncol=8)
  rownames(model) <- c(amplicons, "pooled")
  colnames(model) <- c("mean_cov_target", "mean_cov_reference", "Shapiro_Wilk_W_CN", "Shapiro_Wilk_p_CN", "mean_CN", "sd_CN", "median_CN", "mad_CN")

  model[, "mean_cov_target"] <- c(m.coverage.amplicon[amplicons], mean(coverage.target[amplicons, ]))
  model[, "mean_cov_reference"] <- c(m.coverage.amplicon[amplicons], mean(coverage.source[amplicons, ]))

  model[, "mean_CN"] <- c(apply(CN.source, 1, mean), mean(CN.source))
  model[, "sd_CN"] <- c(apply(CN.source, 1, sd), sd(CN.source))
  model[, "median_CN"] <- c(apply(CN.source, 1, median), median(CN.source))
  s <- c(apply(CN.source, 1, mad), mad(CN.source))
  model[, "mad_CN"] <- s

  for (i in 1:namplicons) {
    cn <- CN.source[i, ]
    shapiro <- shapiro.test(cn)
    model[i, "Shapiro_Wilk_W_CN"] <- shapiro$statistic
    model[i, "Shapiro_Wilk_p_CN"] <- shapiro$p.value
  }


  #2: AMPLICON MATRIX

  CN.a <- CN.target
  P.a <- matrix(nrow=namplicons, ncol=nsamples)
  colnames(P.a) <- samples
  rownames(P.a) <- amplicons


  for (a in amplicons){
    for (s in samples){
      cn <- CN.target[a,s]
      cn <- abs(cn-2)
      if (method.pooled == "pooled") SD <- model["pooled", "mad_CN"]
      if (method.pooled == "amplicon") SD <- model[a, "mad_CN"]
      P.a[a,s] <- pnorm(cn, mean=0, sd=SD, lower.tail=FALSE)
    }
  }

  #3: GENE MATRIX

  genes <- c()
  genescut <- strsplit(rownames(CN.a), "_|-")
  for(i in 1:length(genescut)) {
    genes <- c(genes,genescut[[i]][1])
  }
  genes2 <- unique(genes)
  ngenes <- length(genes2)

  CN.g <- matrix(nrow=ngenes, ncol=nsamples)
  rownames(CN.g) <- genes2
  colnames(CN.g) <- samples
  P.g <- CN.g

  coverages.genes <- rep(NA, ngenes)
  names(coverages.genes) <- genes2

  pfisher <- function(p) {
    q <- -2*sum(log(p))
    df <- 2*length(p)
    p.fisher <- pchisq(q, df, lower.tail=FALSE)
    return(p.fisher)
  }

  genes3 <- c()
  genescut2 <- strsplit(names(scale.amplicons), "_|-")
  for(i in 1:length(genescut2)) {
    genes3 <- c(genes3,genescut2[[i]][1])
  }

  for (g  in genes2) {
    index.g <- which(genes == g)
    CN <- CN.a[index.g, ,drop=FALSE]
    CN.g[g, ] <- apply(CN, 2, mean)
    P <- P.a[index.g, , drop=FALSE]
    P.g[g, ] <- apply(P, 2, pfisher)
    index.g2 <- which(genes3 == g)
    coverages.genes[g] <- paste(paste(names(scale.amplicons[index.g2]), ": ", round(2 * scale.amplicons[index.g2],0), sep=""), collapse=", ")
  }
  attr(CN.target, "scale.sample") <- scale.samples
  attr(CN.g, "coverages.genes") <- coverages.genes
  attr(CN.target, "scale.amplicon") <- scale.amplicons
  CN <- list()

  CN[["model"]] <- model
  CN[["target"]] <- CN.target
  CN[["CN.a"]] <- CN.a
  CN[["P.a"]] <- P.a
  CN[["CN.g"]] <- CN.g
  CN[["P.g"]] <- P.g


  return(CN)
}
