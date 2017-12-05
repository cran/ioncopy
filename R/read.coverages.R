read.coverages <- function(chip.names, file.names, anno.col="Target") {
  n <- length(file.names)
  cat("Reading coverages from", n, "files...\n")
  coverage <- NULL
  namplicons <- 0
  nsamples <- 0 
  if (n > 0) for (i in 1:n) {
    X <- as.matrix(read.table(file.names[i], header=TRUE, sep="\t"))
    m <- ncol(X)
    if (anno.col %in% colnames(X)) {
      rownames(X) <- X[, anno.col]
      if (i == 1) targets <- rownames(X)
      if (length(setdiff(rownames(X), targets)) == 0) {
        index <- NULL
        for (j in 1:m) suppressWarnings(if (as.logical(sum(is.finite(as.numeric(X[, j]))))) index <- c(index, j))
        Y <- apply(X[targets, index], 1:2, as.numeric)
        col.names <- paste(chip.names[i], colnames(Y), sep="_")
        col.names <- sub("_IonXpress", "", sub(".xls", "", sub(".txt", "", col.names)))  
        colnames(Y) <- col.names
        coverage <- cbind(coverage, Y)
        if (i == 1) rownames(coverage) <- targets
      }
      else stop("All samples need to be analyzed using the same sequencing panel")
    }
    else stop(paste("No column ", anno.col, " in ", chip.names[i], sep=""))
    namplicons <- nrow(coverage)
    nsamples <- ncol(coverage)
  }
  cat("Coverage data of", namplicons, "amplicons and", nsamples, "samples\n") 
  return(coverage)
}
