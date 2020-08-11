library(shiny)

shinyUI(fluidPage(

  splitLayout(

    wellPanel(verticalLayout(
      h4(strong("Data Upload")),
      fileInput("target.files", "Target Coverage Files:", multiple=TRUE),
      fileInput("source.files", "Reference Coverage Files:", multiple=TRUE),
      textInput("anno.col", "Annotation Column:", value="Target"),
      textInput("thres.cov", "Minimum Amplicon Coverage:", value="100"),
      hr(),
      h4(strong("Start Analysis")),
      actionButton("run", "Go!")
    )),


    wellPanel(verticalLayout(
      h4(strong("Analysis Parameters")),
      radioButtons("ana.mode",label = "Analysis mode:",choices = c("gene-wise","amplicon-wise"),selected = "gene-wise"),
      checkboxGroupInput("choice", "Multiple Testing Correction:", c("samples" = "samples", "genes/amplicons" = "genes/amplicons"), selected=c("samples", "genes/amplicons"), inline=TRUE),
      radioButtons(inputId="multtest", label = "Method:",choices = c("Bonferroni (FWER)","Benjamini-Hochberg (FDR)", "none"),selected = "Bonferroni (FWER)"),
      textInput("thres.p", label="Significance level:", value="0.05"),
      textInput("sig.call", label="Minimum Number of Significant Amplicons:", value="0"),
      textInput("sig.per", label="Minimum Percentage of Significant Amplicons:", value="0"),
      hr(),
      h4(strong("Heatmap Parameters")),
      textInput("thres.percent", "Minimum percentage:", value="0"),
      checkboxGroupInput("cluster", "Clustering:", c("samples" = "samples", "genes/amplicons" = "genes/amplicons"), selected=c("genes/amplicons","samples" ), inline=TRUE),
      radioButtons(inputId="data", label = "Data Type:",choices = c("copy numbers","CNA calls"),selected = "CNA calls")

    )),


    mainPanel(
      h4(strong("Ioncopy Results")),
      p(
        downloadButton("target.coverage", "Tar. Coverages: matrix"),
        downloadButton("source.coverage", "Ref. Coverages: matrix"),
        downloadButton("model", "QC: list"),
        #downloadButton("model2", "QC_samples: list"),
        downloadButton("CN.a", "CNs Amplicons: matrix"),
        downloadButton("CN.g", "CNs Genes: matrix")
      ),
      p(
        downloadButton("CNA.list", "CNAs: list"),
        downloadButton("GAIN", "Gains: matrix"),
        downloadButton("LOSS", "Losses: matrix"),
        downloadButton("sumCNA.samples", "CNAs samples: list"),
        downloadButton("sumCNA.wise", "CNAs genes/amplicons: list"),
        downloadButton("zip", "All files: zip")
      ),
      p(
      textOutput("error"),
      textOutput("status"),
      textOutput("summary")
      ),
      tags$head(tags$style("#error{color: red; font-size: 16px;}")),
      tags$head(tags$style("#status{color: black; font-size: 16px;}")),
      tags$head(tags$style("#summary{color: black; font-size: 16px;}")),
      plotOutput("heatmap", width=900, height=900, inline=FALSE)
    ),
    cellWidths=c("26%","25%","75%")
  )
))
