
convert <- function(CN, P, type="amplicon") {
  namplicons <- nrow(CN)
  nsamples <- ncol(CN)
  result <- NULL

  for (i in 1:nsamples) {
    res <- cbind(colnames(CN)[i], rownames(CN), CN[,i], P[, i])
    result <- rbind(result, res)
  }
  colnames(result) <- c("sample", type, "CN", "p")
  rownames(result) <- paste(result[, "sample"], result[, type], sep="_")
  return(result)
}



call.CNA <- function(CNA, analysis.mode="gene-wise", method.p="samples_genes/amplicons", method.mt="bonferroni", thres.p=0.05, sig.call=0, sig.per=0) {

    nsamples <- ncol(CNA[["CN.a"]])
    samples <- colnames(CNA[["CN.a"]])

    #analysis.mode:
    if (analysis.mode == "amplicon-wise") A <- convert(CNA[["CN.a"]], CNA[["P.a"]], type = "amplicon")
    if (analysis.mode == "gene-wise") A <- convert(CNA[["CN.g"]], CNA[["P.g"]], type = "gene")


    if (analysis.mode == "amplicon-wise"){
      wise <- "amplicon"
      wiselist <- rownames(CNA[["CN.a"]])
    }
    if (analysis.mode == "gene-wise"){
      wise <- "gene"
      wiselist <- rownames(CNA[["CN.g"]])
    }



    if (method.p == "") name.p <- "p_raw"
    else name.p <- paste("p", method.p, method.mt, sep="_")


    if (method.p == "") name.call <- paste("call", thres.p, sep="_")
    else name.call <- paste("call", method.p, method.mt, thres.p, sep="_")

    A <- cbind(A, NA, "NORMAL")

    colnames(A)[5] <- name.p


    #for only amplicons/genes:
    if (method.p == "genes/amplicons") {

      if(method.mt == "none") {
        A[,name.p] <- A[,"p"]
      }
      else{
      for (x in wiselist) {
        index <- which(A[, wise] == x)
        p <- as.numeric(A[index, "p"])
        A[index, name.p] <- p.adjust(p, method=method.mt)
      }
      }
    }

    #for only samples:
    if (method.p == "samples"){

      if(method.mt == "none") {
        A[,name.p] <- A[,"p"]
      }
      else{
      for (y in samples) {
        index <- which(A[, "sample"] == y)
        p <- as.numeric(A[index, "p"])
        A[index, name.p] <- p.adjust(p, method=method.mt)
      }
      }
    }

    #for both:
    if (method.p == "samples_genes/amplicons") {

      if(method.mt == "none") {
        A[,name.p] <- A[,"p"]
      }
      else{
          p <- as.numeric(A[, "p"])
          A[, name.p] <- p.adjust(p, method=method.mt)
      }

    }

    if(method.p == "") {
      A[,name.p] <- A[,"p"]
    }

    #calls: thres.p
    colnames(A)[6] <- name.call
    index.p <- which(as.numeric(A[, name.p]) < thres.p)
    index.gain <- which(as.numeric(A[, "CN"]) > 2)
    index.gain <- intersect(index.gain, index.p)
    index.loss <- which(as.numeric(A[, "CN"]) < 2)
    index.loss <- intersect(index.loss, index.p)
    if (length(index.gain) > 0) A[index.gain, name.call] <- "GAIN"
    if (length(index.loss) > 0) A[index.loss, name.call] <- "LOSS"


    if (analysis.mode == "gene-wise"){
      A <- cbind(A," ")
      colnames(A)[7] <- "ncalls"
      amplicons <- rownames(CNA[["CN.a"]])
      CN.a <- CNA[["CN.a"]]
      P.a <- CNA[["P.a"]]
      liste <- which(A[,name.call]!="NORMAL")
          for(l in liste){
            w <- A[l,"gene"]
            s <- A[l, "sample"]
            d <- A[l, name.call]
            a <- amplicons[grep(w,amplicons)]
            cn <- CN.a[a,s]
            p <- P.a[a,s]
            if (d == "GAIN") index.cn <- which(cn > 2)
            else index.cn <- which(cn < 2)
            index.p <- which(p < 0.05)
            index <- intersect(index.cn, index.p)
            result <- paste(length(index), length(cn), sep = " of ")
            if((length(index) >= sig.call)&&((length(index)/length(a))*100>=sig.per)){
                A[l,"ncalls"] <- result
            }
            else A[l, name.call] <- "NORMAL"
            }
    }




    if (analysis.mode == "amplicon-wise") nwise <- nrow(CNA[["CN.a"]])
    if (analysis.mode == "gene-wise") nwise <- nrow(CNA[["CN.g"]])

    LOSS <- matrix(nrow=nwise, ncol=nsamples)
    LOSS <- apply(LOSS, c(1, 2), function(x) 0)
    GAIN <- matrix(nrow=nwise, ncol=nsamples)
    GAIN <- apply(GAIN, c(1, 2), function(x) 0)

    if (analysis.mode == "amplicon-wise"){
      rownames(LOSS) <- rownames(CNA[["CN.a"]])
      colnames(LOSS) <- colnames(CNA[["CN.a"]])
      rownames(GAIN) <- rownames(CNA[["CN.a"]])
      colnames(GAIN) <- colnames(CNA[["CN.a"]])
    }
    if (analysis.mode == "gene-wise"){
      rownames(LOSS) <- rownames(CNA[["CN.g"]])
      colnames(LOSS) <- colnames(CNA[["CN.g"]])
      rownames(GAIN) <- rownames(CNA[["CN.g"]])
      colnames(GAIN) <- colnames(CNA[["CN.g"]])
    }



     for (i in 1:nsamples){

     index <- ((i-1)*nwise+1):(i*nwise)

     index.g <- which(A[index, name.call]=="GAIN")
     index.l <- which(A[index, name.call]=="LOSS")
     GAIN[index.g, i] <- 1
     LOSS[index.l, i] <- 1
     }


  calls <- list()
  calls[["tab"]] <- A
  calls[["loss"]] <- LOSS
  calls[["gain"]] <- GAIN
  calls[["CN.a"]] <- CNA[["target"]]

  if (analysis.mode == "amplicon-wise") calls[["CN"]] <- CNA[["CN.a"]]
  if (analysis.mode == "gene-wise") calls[["CN"]] <- CNA[["CN.g"]]


  return(calls)

  }


