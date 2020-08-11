library(shiny)
library(zip)

options(shiny.maxRequestSize=100*1024^2)

write.xls <- function(x, file, colnames.1="") {
  y <- cbind(rownames(x), x)
  colnames(y)[1] <- colnames.1
  write.table(y, file, row.names=FALSE, col.names=TRUE, sep="\t")
}

shinyServer(function(input, output) {

  target.files <- reactive(input$target.files)
  source.files <- reactive(input$source.files)
  anno.col <-  reactive(input$anno.col)
  thres.cov <- reactive(as.numeric(input$thres.cov))

  thres.p.cna <- reactive(as.numeric(input$thres.p))
  method.p.cna <- reactive(paste(input$choice, collapse="_"))
  analysis.mode.cna <- reactive(input$ana.mode)
  method.mt.cna <- reactive(input$multtest)
  sig.amp <- reactive(as.numeric(input$sig.call))
  sig.per <- reactive(as.numeric(input$sig.per))

  thres.percent <-  reactive(as.numeric(input$thres.percent))
  cluster.genes <- reactive("genes/amplicons" %in% input$cluster)
  cluster.samples <- reactive("samples" %in% input$cluster)
  heatdata <- reactive(input$data)

  observeEvent(input$run, {withProgress(message="", value=0, {

    output$heatmap <- renderPlot(NULL)
    output$status <- renderText(NULL)

    incProgress(1/5, message="Reading files...")
    file.names <- target.files()[, "datapath"]
    chip.names <- substr(tools::file_path_sans_ext(target.files()[, "name"]), 1, 40)
    coverage.target <- try(read.coverages(chip.names, file.names, anno.col()))
    #print(colnames(coverage.target))
    file.names <- source.files()[, "datapath"]
    chip.names <- substr(tools::file_path_sans_ext(source.files()[, "name"]), 1, 40)
    coverage.source <- try(read.coverages(chip.names, file.names, anno.col()))
    if (is.null(coverage.source)) {
      coverage.source <- coverage.target
      one.cohort <- TRUE
    }
    else one.cohort <- FALSE

    ngenes <- length(unique(sapply(rownames(coverage.target), function(x){strsplit(x, "_|-")[[1]][1]})))

    output$error <- renderText({
      validate(
        need(!is.null(coverage.target), "ERROR: Please specify files for target cohort!"),
        need(!class(coverage.target) == "try-error", paste("ERROR (target cohort): ", attr(coverage.target, "condition")$message, sep="")),
        need((!class(coverage.source) == "try-error" || one.cohort), paste("ERROR (reference cohort): ", attr(coverage.source, "condition")$message, sep="")),
        need(length(setdiff(rownames(coverage.target), rownames(coverage.source))) == 0, "ERROR: Reference and target samples need to be analyzed using the same sequencing panel")
      )
      output$status <- renderText(paste(nrow(coverage.target), "amplicons, ", ngenes, "genes, ", ncol(coverage.target), "samples (target cohort), ", ncol(coverage.source), "samples (reference cohort)."))
      return(NULL)
    })


    #Die Vor-Matritzen:
    try ({
    coverage.source <- coverage.source[rownames(coverage.target), ]

    coverage.source2 <- t(coverage.source)
    coverage.target2 <- t(coverage.target)

    output$target.coverage <- downloadHandler(
      filename = function(){"coverage_target.xls"},
      content = function(file){write.xls(coverage.target2, file, colnames.1="amplicon")}
    )
    output$source.coverage <- downloadHandler(
      filename = function(){"coverage_reference.xls"},
      content = function(file){write.xls(coverage.source2, file, colnames.1="amplicon")}
    )

    incProgress(1/5, message="Fitting model...")
    cn <- assess.CNA(coverage.target, coverage.source, method.pooled="amplicon", thres.cov=thres.cov())

    CN.a <- t(cn[["CN.a"]])
    CN.g <- t(cn[["CN.g"]])

    output$model <- downloadHandler(
      filename = function(){"QC_amplicons.xls"},
      content = function(file){write.xls(cn[["model"]], file, colnames.1="amplicon")}
    )
    #output$model2 <- downloadHandler(
      #filename = function(){"QC_samples.xls"},
      #content = function(file){write.xls(cn[["model2"]], file, colnames.1="samples")}
    #)
    output$CN.a <- downloadHandler(
      filename = function(){"CN_amplicons.xls"},
      content = function(file){write.xls(CN.a, file, colnames.1="amplicons")}
    )
    output$CN.g <- downloadHandler(
      filename = function(){"CN_genes.xls"},
      content = function(file){write.xls(CN.g, file, colnames.1="genes")}
    )


    if(analysis.mode.cna()=="amplicon-wise"){
      wise <- "amplicon"
      variable <- "amplicons"
    }
    if(analysis.mode.cna()=="gene-wise") {
      wise <- "gene"
      variable <- "genes"
    }

    incProgress(1/5, message="Calling amplifications/genes...")

    if (method.mt.cna() == "Bonferroni (FWER)") method.mt <- "bonferroni"
    if (method.mt.cna() == "Benjamini-Hochberg (FDR)") method.mt <- "BH"
    if (method.mt.cna() == "none") {method.mt <- "none"}
    cna <- call.CNA(cn, analysis.mode=analysis.mode.cna(), method.p=method.p.cna(), method.mt=method.mt, thres.p=thres.p.cna(), sig.call=sig.amp(), sig.per=sig.per())


    sumcna <- summarize.CNA(cna)

    Loss <- cna[["loss"]]
    Gain <- cna[["gain"]]
    ngains <- length(which(Gain == "1"))
    nloss <- length(which(Loss == "1"))
    nwise <- nrow(Loss)
    nsamples <- ncol(Loss)
    dim <- nsamples*nwise
    percg <- paste(round((ngains/dim)*100,1), "%", sep = "")
    percl <- paste(round((nloss/dim)*100,1), "%", sep = "")

    output$summary <- renderText(paste("Detection of ", ngains," (", percg, ") gains and ", nloss, " (", percl, ") losses in ", dim, " measurements.", sep = ""))


    output$CNA.list <- downloadHandler(
      filename = function(){"CNA.xls"},
      content = function(file){write.xls(cna[["tab"]], file,colnames.1=paste("sample",wise, sep = "_"))}
    )

    Losses <- t(Loss)
    Gains <- t(Gain)

    output$LOSS <- downloadHandler(
      filename = function(){"LOSSES.xls"},
      content = function(file){write.xls(Losses, file, colnames.1=wise)}
    )

    output$GAIN <- downloadHandler(
      filename = function(){"GAINS.xls"},
      content = function(file){write.xls(Gains, file, colnames.1=wise)}
    )

    output$sumCNA.samples <- downloadHandler(
      filename = function(){"CNA_samples.xls"},
      content = function(file){write.xls(sumcna[["samples"]], file, colnames.1="samples")}
    )

    output$sumCNA.wise <- downloadHandler(
      filename = function(){paste("CNA_",variable,".xls", sep = "")},
      content = function(file){write.xls(sumcna[[wise]], file, colnames.1=wise)}
    )

    output$heatmap <- renderPlot(
      heatmap.CNA(cna, thres.percent=thres.percent(), cluster.genes=cluster.genes(), cluster.samples=cluster.samples(), type=heatdata()),
      width=900,
      height=900,
      res=200,
      pointsize=24
    )

    # All files in a zip folder:
    # Second function for output heatmaps:
    heatmap.CNA2 <- function(CNA, fileName, thres.percent=1, cluster.genes=TRUE, cluster.samples=TRUE, type="CNA calls" ,method.dist="manhattan", method.link="average", mar=3, cex=0.7) {



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
            pos.legend <- 1.14
          }
          else { pos.legend <- 1.18 }

          pdf(file = fileName, width=400, height=400, paper = "a4r", pointsize = 36)
          heatmap(GL, Rowv=Rowv, Colv=Colv, scale="none", distfun=function(x){dist(x, method=method.dist)}, hclustfun=function(x){hclust(x, method=method.link)}, reorderfun=function(d, w){reorder(d, w, agglo.FUN=sum)}, breaks=c(-1.5, -0.5, 0.5, 1.5, 10), col=c("green", "black", "red", "blue"), margins=c(mar, mar), cexRow=cexRow, cexCol=cexCol,
                  legend(ncol(GL)*pos.legend, 0.5, yjust = 0, legend = c("gain", "loss", "normal"), col = c("red", "green", "black"), pch = 15, cex=0.3, xpd = NA))
          dev.off()
         }

        #return(GL)
      }



      #Copy Numbers

      if ( type == "copy numbers"){

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


        CN <- CNA[["CN"]]
        CN <- CN[index,]

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
        pdf(file = fileName, width=400, height=400, paper = "a4r", pointsize = 36)
        heatmap(CN, Rowv=Rowv, Colv=Colv, scale="none", distfun=function(x){dist(x, method=method.dist)}, hclustfun=function(x){hclust(x, method=method.link)}, reorderfun=function(d, w){reorder(d, w, agglo.FUN=sum)}, breaks=c(0:5,10,10^5), col=c("darkgreen", rgb(0, .15, 0), rgb(.15, 0, 0), "darkred", "red", "orange", "yellow"), margins=c(mar, mar), cexRow=cexRow, cexCol=cexCol,
                add.expr= legend(ncol(CN)*pos.legend, 0.5, yjust = 0, legend = c("CN < 1", "CN < 2", "CN > 2", "CN > 3", "CN > 4", "CN > 5", "CN > 10"), col = c("darkgreen", rgb(0, .15, 0), rgb(.15, 0, 0), "darkred", "red", "orange", "yellow"), pch = 15, cex = 0.3, xpd=NA))
        dev.off()
        #return(CN)
      }


    }
    chip.names <- substr(tools::file_path_sans_ext(target.files()[, "name"]), 1, 40)
    output$zip <- downloadHandler(
      filename = function(){
        paste0(paste("ioncopy", chip.names, sep = "_"),".zip")

      },
      content = function(file){
        #go to a temp dir to avoid permission issues
        tmpdir <- tempdir()
        setwd(tempdir())
        files <- c()

        fileName <- "coverage_target.xls"
        files <- c(fileName,files)
        write.xls(coverage.target2, fileName, colnames.1="amplicon")
        fileName <- "coverage_reference.xls"
        files <- c(fileName,files)
        write.xls(coverage.source2, fileName, colnames.1="amplicon")
        fileName <- "QC_amplicons.xls"
        files <- c(fileName,files)
        write.xls(cn[["model"]], fileName, colnames.1="amplicon")
        fileName <- "CN_amplicons.xls"
        files <- c(fileName,files)
        write.xls(CN.a, fileName, colnames.1="amplicons")
        fileName <- "CN_genes.xls"
        files <- c(fileName,files)
        write.xls(CN.g, fileName, colnames.1="genes")
        fileName <- "CNA.xls"
        files <- c(fileName,files)
        write.xls(cna[["tab"]], fileName,colnames.1=paste("sample",wise, sep = "_"))
        fileName <- "LOSSES.xls"
        files <- c(fileName,files)
        write.xls(Losses, fileName, colnames.1=wise)
        fileName <- "GAINS.xls"
        files <- c(fileName,files)
        write.xls(Gains, fileName, colnames.1=wise)
        fileName <- "CNA_samples.xls"
        files <- c(fileName,files)
        write.xls(sumcna[["samples"]], fileName, colnames.1="samples")
        fileName <- paste("CNA_",variable,".xls", sep = "")
        files <- c(fileName,files)
        write.xls(sumcna[[wise]], fileName, colnames.1=wise)
        fileName <- "Heatmap_copy_numbers.pdf"
        files <- c(fileName,files)
        heatmap.CNA2(cna, fileName, thres.percent=thres.percent(), cluster.genes=cluster.genes(), cluster.samples=cluster.samples(), "copy numbers")
        fileName <- "Heatmap_CNA_calls.pdf"
        files <- c(fileName,files)
        heatmap.CNA2(cna, fileName, thres.percent=thres.percent(), cluster.genes=cluster.genes(), cluster.samples=cluster.samples(), type="CNA calls")

        #create the zip file
        zipr(file,files)
      },
      contentType = "application/zip"
    )


  })
  })
  })
})
