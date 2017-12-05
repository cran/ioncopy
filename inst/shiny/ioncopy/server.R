library(shiny)

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
    output$error <- renderText({
      validate(
        need(!is.null(coverage.target), "ERROR: Please specify files for target cohort!"),
        need(!class(coverage.target) == "try-error", paste("ERROR (target cohort): ", attr(coverage.target, "condition")$message, sep="")),
        need((!class(coverage.source) == "try-error" || one.cohort), paste("ERROR (reference cohort): ", attr(coverage.source, "condition")$message, sep="")),
        need(length(setdiff(rownames(coverage.target), rownames(coverage.source))) == 0, "ERROR: Reference and target samples need to be analyzed using the same sequencing panel")
      )
      output$status <- renderText(paste(nrow(coverage.target), "amplicons, ", ncol(coverage.target), "samples (target cohort), ", ncol(coverage.source), "samples (reference cohort)."))
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
    cn <- assess.CNA(coverage.target, coverage.source, method.pooled="amplicon")

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

    if (method.mt.cna() == "Bonferroni (FWER)") method.mt <- "Bonferroni"
    if (method.mt.cna() == "Benjamini-Hochberg (FDR)") method.mt <- "BH"
    if (method.mt.cna() == "none") {method.mt <- "none"}
    cna <- call.CNA(cn, analysis.mode=analysis.mode.cna(), method.p=method.p.cna(), method.mt=method.mt, thres.p=thres.p.cna(), sig.call=sig.amp())


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
  })
  })
  })
})
