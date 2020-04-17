### -------------------------------------------------------------------- ###
###                            MODULE DIFFEXP
### -------------------------------------------------------------------- ###
library(edgeR)

### FUNCTIONS
### -------------------------------------------------------------------- ###


### UI
### -------------------------------------------------------------------- ###
computeDiffExpOutput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  tagList(
    fluidRow(column(12, verbatimTextOutput(ns("deSummary")))),
    br(),
    fluidRow(
      column(6, plotOutput(ns("dePval"))),
      column(6, plotOutput(ns("deFdr")))
    ),
    br(),
    # fluidRow(
    #   column(6, plotOutput(ns("deSmear")))
      # column(6, plotOutput(ns("deFdr")))
    # )
  )
}


### SERVER
### -------------------------------------------------------------------- ###
computeDiffExp <- function(input, output, session, button, dgeObjNorm, feature, condition) {

  # Reactive values
  rv <- reactiveValues(
    dgeDecide = NULL,
    lrtModel = NULL,
    dgeDf = NULL,
    dgeChangedGenes = NULL
  )
  
  # Compute the differentially expressed genes
  observeEvent(button(), {
    # Compute the differential expression of genes, for a control vs treatment design
    ## prepare dge object for differential expression (DE) computation
    dgeObjNorm <- estimateCommonDisp(dgeObjNorm())
    dgeObjNorm <- estimateGLMTrendedDisp(dgeObjNorm)
    dgeObjNorm <- estimateTagwiseDisp(dgeObjNorm)

    ## get design matrix
    features <- dgeObjNorm$samples[, as.character(feature())]
    features <- features == as.character(condition())
    design <- model.matrix(~ features)

    ## deduce LRT model and deduce differentially expressed genes
    fit <- glmFit(dgeObjNorm, design)
    rv$lrtModel <- glmLRT(fit, coef = 2)
    rv$dgeDf <- as.data.frame(topTags(rv$lrtModel, n = Inf))
    rv$dgeDecide <- decideTestsDGE(rv$lrtModel, p = 0.05, adjust = "BH")
    rv$dgeChangedGenes <- row.names(dgeObjNorm()[as.logical(rv$dgeDecide), ])
  })
  
  # Plot the the results
  ## plot summary
  output$deSummary <- renderPrint({
    if(!is.null(rv$dgeDecide)){
      summary(rv$dgeDecide)
      }
  })
  ## plot p-values histogram
  output$dePval <- renderPlot({
    if(!is.null(rv$dgeDecide)){
      hist(rv$dgeDf$PValue,
           las = 2,
           main = "Histogram of p-values",
           xlab = "p-values")
    } else{}
  })
  ## plot adjusted p-values histogram
  output$deFdr <- renderPlot({
    if(!is.null(rv$dgeDf)){
      hist(rv$dgeDf$FDR,
           las = 2,
           main = "Histogram of FDR",
           xlab = "FDR")
      abline(v = 0.05, col = "blue")
    }
  })

  # Return results
  return(rv)
}