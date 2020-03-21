### -------------------------------------------------------------------- ###
###                            MODULE DGEOBJ
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
    fluidRow(verbatimTextOutput(ns("deSummary"))),
    br(),
    fluidRow(
      column(6, plotOutput(ns("dePval"))),
      column(6, plotOutput(ns("deFdr")))
    )
  )
}


### SERVER
### -------------------------------------------------------------------- ###
computeDiffExp <- function(input, output, session, button, dgeObjNorm, feature, condition) {

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

    ## deduce LRT model
    fit <- glmFit(dgeObjNorm, design)
    lrtModel <- glmLRT(fit, coef = 2)
    dgeDf <- as.data.frame(topTags(lrtModel, n = Inf))
    dgeDecide <- decideTestsDGE(lrtModel, p = 0.05, adjust = "BH")
    print("l43")
    
    # Plot the the results
    ## plot summary
    output$deSummary <- renderPrint({
      summary(dgeDecide)
    })
    ## plot the pvalues and FDR histograms
    output$dePval <- renderPlot({
      hist(dgeDf$PValue,
           las = 2,
           main = "Histogram of p-values",
           xlab = "p-values")
    })
    output$deFdr <- renderPlot({
      hist(dgeDf$FDR,
           las = 2,
           main = "Histogram of FDR",
           xlab = "FDR")
      abline(v = 0.05, col = "blue")
    })
    
    # Return results
    results <- list(lrtModel, dgeDf, dgeDecide)
    names(results) <- c("lrtModel", "dgeDf", "dgeDecide")
    return(results)
  })
}