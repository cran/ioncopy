
qc.CN <- function(CN, methods=c("sd", "mad"), variables=c("samples", "genes")) {
  nmethods <- length(methods)
  samples <- colnames(CN)
  nsamples <- length(samples)
  genes <- sapply(strsplit(rownames(CN), "_|-"), function(x){x[1]})
  genes.u <- sort(unique(genes))
  ngenes <- length(genes.u)
  n.u <- sapply(genes.u, function(g){length(which(genes == g))})
  index <- which(n.u >= 2)

  if (variables[1] == "samples") {
    QC <- matrix(nrow=nmethods, ncol=nsamples)
    rownames(QC) <- methods
    colnames(QC) <- samples
    variables.i <- 2
  }

  if (length(index) > 0) {
    genes.qc <- genes.u[index]
    n <- length(genes.qc)

    if (variables[1] == "genes") {
      QC <- matrix(nrow=nmethods, ncol=length(genes.qc))
      rownames(QC) <- methods
      colnames(QC) <- genes.qc
      variables.i <- 1
    }

    Q <- list()
    for (m in methods) {
      Q[[m]] <- matrix(nrow=n, ncol=nsamples)
      rownames(Q[[m]]) <- genes.qc
      colnames(Q[[m]]) <- samples
      for (g in genes.qc) {
        index <- which(genes == g)
        cn <- CN[index, ]
        #print(cn)
        if (m == "sd") qc <- apply(cn, 2, sd)
        if (m == "mad") qc <- apply(cn, 2, mad)
        Q[[m]][g, ] <- qc
      }
    }
    for (m in methods) {
      if (m == "sd") QC[m, ] <- apply(Q[[m]], variables.i, mean)
      if (m == "mad") QC[m, ] <- apply(Q[[m]], variables.i, median)
    }
  }
  return(QC)
}


summarize.CNA <- function(calls){

  samples <- colnames(calls[["CN"]])
  nsamples <- ncol(calls[["CN"]])
  wise <- rownames(calls[["CN"]])
  nwise <- nrow(calls[["CN"]])

  Gain <- calls[["gain"]]
  Loss <- calls[["loss"]]
  ngains <- length(which(Gain=="1"))
  nlosses <- length(which(Loss == "1"))


  X1 <- calls[["tab"]]
  out <- which(X1[,6]=="NORMAL")
  X <- X1[-out,]
  wisename <- colnames(X)[2]

  cna.samples <- matrix(nrow=nsamples, ncol=4)
  rownames(cna.samples) <- samples
  colnames(cna.samples)<- c("coverage_median","intra_gene_inconsistency", "n_alterations", "percent_alterations")
  nalt <- sapply(samples, function(s){ length(which(s==X[,"sample"]))})
  palt <- paste(round((nalt/nwise)*100,1), sep="")
  cna.samples[,"n_alterations"] <- nalt
  cna.samples[, "percent_alterations"] <- palt
  coverages.samples <- attr(calls[["CN.a"]], "scale.sample")
  for(i in samples){
    cna.samples[i, "coverage_median"] <- round(coverages.samples[i],0)
  }

  QC <- qc.CN(calls[["CN.a"]], methods = c("sd"), variables = c("samples"))
  for(i in samples){
    cna.samples[i, "intra_gene_inconsistency"] <- round(QC[,i],2)
  }

    cna.samples <- cbind(cna.samples, "", "", "", "")
    colnames(cna.samples)[(length(colnames(cna.samples))-3):length(colnames(cna.samples))] <- c("n_gains", "percent_gains", "n_losses", "percent_losses")
    for(i in samples){
      numberl <- length(which(Loss[, i]=="1"))
      cna.samples[i,"n_losses"] <- numberl
      numberg <- length(which(Gain[, i]=="1"))
      cna.samples[i,"n_gains"] <- numberg
      cna.samples[i, "percent_gains"] <- paste(round((numberg/ngains*100),1), sep = "")
      cna.samples[i, "percent_losses"] <- paste(round((numberl/nlosses*100),1), sep = "")
    }

    cna.samples <- cbind(cna.samples, "")
    colnames(cna.samples)[9] <- "alterations"
    for(s in samples){
      index <- which(s == X[, "sample"])
      w <- X[index, wisename]
      cn <- round(as.numeric(X[index, "CN"]), 1)
      p <- signif(as.numeric(X[index, "p"]), 2)
      ind <- order(w)
      w <- w[ind]
      cn <- cn[ind]
      p <- p[ind]

    if(wisename=="gene"){
      dec2 <- c()
      dec <- X[index, 7]
      dec <- dec[ind]
      for(d in dec){
        split <- strsplit(d, " of ")
        dec2 <- c(dec2, paste(split[[1]][1], split[[1]][2], sep = "/"))
      }
      dec <- dec2
      gl <- X[index, 6]
      gl <- gl[ind]
      x <- paste(w," (", gl, ", CN=", cn, ", ", "p=", p,", ", "n_detected=", dec,")", sep = "")
      cna.samples[s,"alterations"] <- paste(x, collapse = "; ")
      if(cna.samples[s, "n_alterations"] == "0") cna.samples[s,"alterations"] <- paste("")
    }

    if(wisename =="amplicon"){
      gl <- X[index, 6]
      gl <- gl[ind]
      x <- paste(w, " (", gl, ", CN=", cn, ", ", "p=", p,")", sep = "")
      cna.samples[s,"alterations"] <- paste(x, collapse = "; ")
      if(cna.samples[s, "n_alterations"] == "0") cna.samples[s,"alterations"] <- paste("")
    }
  }

  for(i in ncol(cna.samples[,"alterations"])){
    sort(cna.samples[i, "alterations"])
  }

  cna.wise <- matrix(nrow=nwise, ncol=2)
  rownames(cna.wise) <- wise
  colnames(cna.wise)<- c("n_alterations", "percent_alterations")
  nalt <- sapply(wise, function(w){ length(which(w==X[,2]))})
  palt <- paste(round((nalt/nsamples)*100,1), sep = "")
  cna.wise[,"n_alterations"] <- nalt
  cna.wise[, "percent_alterations"] <- palt

  if(wisename=="gene") {
    cna.wise <- cbind(0, "",cna.wise)
    colnames(cna.wise)[1] <- "coverage_median"
    colnames(cna.wise)[2] <- "intra_gene_inconsistency"
    coverages.genes <- attr(calls[["CN"]], "coverages.genes")
    for(i in wise){
      cna.wise[i, "coverage_median"] <- coverages.genes[i]
    }
    qc <- qc.CN(calls[["CN.a"]], methods = c("sd"), variables = c("genes"))
    for(i in colnames(qc)){
      cna.wise[i, "intra_gene_inconsistency"] <- round(qc[,i],2)
    }
  }

  if(wisename=="amplicon"){
    cna.wise <- cbind("", cna.wise)
    colnames(cna.wise)[1] <- "coverage_median"
    coverages.amplicons <- attr(calls[["CN.a"]], "scale.amplicon")
    for(i in wise){
      cna.wise[i, "coverage_median"] <- round(as.numeric(2*coverages.amplicons[i]),0)
    }
  }

  cna.wise <- cbind(cna.wise, "", "", "", "")
  colnames(cna.wise)[(length(colnames(cna.wise))-3):length(colnames(cna.wise))] <- c("n_gains", "percent_gains", "n_losses", "percent_losses")
  for(i in wise){
    numberl <- length(which(Loss[i,]=="1"))
    cna.wise[i,"n_losses"] <- numberl
    numberg <- length(which(Gain[i,]=="1"))
    cna.wise[i,"n_gains"] <- numberg
    cna.wise[i, "percent_gains"] <- paste(round((numberg/ngains*100),1), sep = "")
    cna.wise[i, "percent_losses"] <- paste(round((numberl/nlosses*100),1), sep = "")
  }

  cna.wise <- cbind(cna.wise, "")
  colnames(cna.wise)[length(colnames(cna.wise))] <- "alterations"

  for(w in wise){
    index <- which(w == X[, wisename])
    sample <- X[index, "sample"]
    cn <- round(as.numeric(X[index, "CN"]), 1)
    p <- signif(as.numeric(X[index, "p"]), 2)
    ind <- order(sample)
    sample <- sample[ind]
    cn <- cn[ind]
    p <- p[ind]

    if(wisename=="gene") {
      dec2 <- c()
      dec <- X[index, 7]
      dec <- dec[ind]
      for(d in dec){
        split <- strsplit(d, " of ")
        dec2 <- c(dec2, paste(split[[1]][1], split[[1]][2], sep = "/"))
      }
      dec <- dec2
      gl <- X[index, 6]
      gl <- gl[ind]
      x <- paste(sample, " (", gl, ", CN=", cn, ", ", "p=", p,", ", "n_detected=", dec, ")", sep = "")
      cna.wise[w,"alterations"] <- paste(x, collapse = "; ")
      if(cna.wise[w, "n_alterations"] == "0") cna.wise[w,"alterations"] <- paste("")
    }

    if(wisename =="amplicon") {
      gl <- X[index, 6]
      gl <- gl[ind]
      x <- paste(sample, " (", gl, ", CN=", cn, ", ", "p=", p,")", sep = "")
      cna.wise[w,"alterations"] <- paste(x, collapse = "; ")
      if(cna.wise[w, "n_alterations"] == "0") cna.wise[w,"alterations"] <- paste("")
    }

  }

  for(i in ncol(cna.wise[,"alterations"])){
    sort(cna.wise[i, "alterations"])
  }


  sumcna <- list()
  sumcna[["samples"]] <- cna.samples
  if(wisename =="gene") sumcna[["gene"]] <- cna.wise
  if(wisename =="amplicon") sumcna[["amplicon"]] <- cna.wise

  return(sumcna)
}
