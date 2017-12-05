heatmap.CNA <- function(CNA, thres.percent=1, cluster.genes=TRUE, cluster.samples=TRUE, type="CNA calls" ,method.dist="manhattan", method.link="average", mar=3, cex=0.50) {



  #CNA calls

  if ( type == "CNA calls"){


  Loss <- CNA[["loss"]]
  Gain <- CNA[["gain"]]

  nwise <- nrow(Loss)
  nsamples <- ncol(Loss)


  GL <- matrix(nrow=nwise, ncol=nsamples)
  rownames(GL) <- rownames(Loss)
  colnames(GL) <- colnames(Loss)

  GL <- Gain - Loss

  m <- apply(abs(GL), 1, sum)
  index <- which(m/ncol(GL)*100 >= thres.percent)
  GL <- GL[index, ]

  if (length(index) >= 2) {
    if (cluster.genes) {
      Rowv <- as.dendrogram(hclust(dist(GL, method=method.dist), method=method.link))
      Rowv <- reorder(Rowv, -rowMeans(GL[order.dendrogram(Rowv), ]))
    }
    else Rowv <- NA
    if (cluster.samples) {
      Colv <- as.dendrogram(hclust(dist(t(GL), method=method.dist), method=method.link))
      Colv <- reorder(Colv, colMeans(GL[, order.dendrogram(Colv)]))
    }
    else Colv <- NA

    cexRow <- cex * 15 / nrow(GL)
    if (cexRow > cex/2) cexRow <- cex/2
    if (cexRow < 0.1) cexRow <- 0.1
    cexCol <- cex * 15 / ncol(GL)
    if (cexCol > cex/2) cexCol <- cex/2
    if (cexCol < 0.1) cexCol <- 0.1

    if((cluster.genes==FALSE) && (cluster.samples==FALSE)){
      pos.legend <- 1.15
    }
    else { pos.legend <- 1.2 }

    heatmap(GL, Rowv=Rowv, Colv=Colv, scale="none", distfun=function(x){dist(x, method=method.dist)}, hclustfun=function(x){hclust(x, method=method.link)}, reorderfun=function(d, w){reorder(d, w, agglo.FUN=sum)}, breaks=c(-1.5, -0.5, 0.5, 1.5, 10), col=c("green", "black", "red", "blue"), margins=c(mar, mar), cexRow=cexRow, cexCol=cexCol,
            legend(ncol(GL)*pos.legend, 0.5, yjust = 0, legend = c("gain", "loss", "normal"), col = c("red", "green", "black"), pch = 15, cex=0.25, xpd = NA))
  }

  return(GL)
  }



  #Copy Numbers

  if ( type == "copy numbers"){

    Loss <- CNA[["loss"]]
    Gain <- CNA[["gain"]]

    nwise <- nrow(Loss)
    nsamples <- ncol(Loss)

      CN <- CNA[["CN"]]

      if (cluster.genes) {
        Rowv <- as.dendrogram(hclust(dist(CN, method=method.dist), method=method.link))
        Rowv <- reorder(Rowv, -rowMeans(CN[order.dendrogram(Rowv), ]))
      }
      else Rowv <- NA
      if (cluster.samples) {
        Colv <- as.dendrogram(hclust(dist(t(CN), method=method.dist), method=method.link))
        Colv <- reorder(Colv, colMeans(CN[, order.dendrogram(Colv)]))
      }
      else Colv <- NA

      cexRow <- cex * 15 / nrow(CN)
      if (cexRow > cex/2) cexRow <- cex/2
      if (cexRow < 0.1) cexRow <- 0.1
      cexCol <- cex * 15 / ncol(CN)
      if (cexCol > cex/2) cexCol <- cex/2
      if (cexCol < 0.1) cexCol <- 0.1

      if((cluster.genes==FALSE) && (cluster.samples==FALSE)){
        pos.legend <- 1.15
      }
      else { pos.legend <- 1.2 }

      heatmap(CN, Rowv=Rowv, Colv=Colv, scale="none", distfun=function(x){dist(x, method=method.dist)}, hclustfun=function(x){hclust(x, method=method.link)}, reorderfun=function(d, w){reorder(d, w, agglo.FUN=sum)}, breaks=c(0:5,10,10^5), col=c("darkgreen", rgb(0, .15, 0), rgb(.15, 0, 0), "darkred", "red", "orange", "yellow"), margins=c(mar, mar), cexRow=cexRow, cexCol=cexCol,
      add.expr= legend(ncol(CN)*pos.legend, 0.5, yjust = 0, legend = c("CN < 1", "CN < 2", "CN > 2", "CN > 3", "CN > 4", "CN > 5", "CN > 10"), col = c("darkgreen", rgb(0, .15, 0), rgb(.15, 0, 0), "darkred", "red", "orange", "yellow"), pch = 15, cex = 0.25, xpd=NA))
     return(CN)
  }


}


