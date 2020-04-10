### -------------------------------------------------------------------- ###
###                            MODULE DGEOBJ
### -------------------------------------------------------------------- ###
library(edgeR)


### FUNCTIONS
### -------------------------------------------------------------------- ###


### UI
### -------------------------------------------------------------------- ###
computeDgeObjOutput <- function(id) {
  
  ns <- NS(id)
    
  tagList(
    fluidRow(column(6, plotOutput(ns("libSizes"))),
             column(6, plotOutput(ns("logCpmNormDistrib")))
    ),
  )
}


### SERVER
### -------------------------------------------------------------------- ###
computeDgeObj <- function(input, output, session, button, counts, features, annot) {
  
  # Create dge object for RNA-Seq analysis
  ## create dge object (an EdgeR object) from a RNA-seq matrix
  dgeObj <- eventReactive(button(), {
    if(is.null(annot)){
      DGEList(counts[-1,], samples = features)
    } else {
      DGEList(counts[-1,], samples = features, genes = annot[,2])
    }
  })
  
  ## normalize counts  
  dgeObjNorm <- eventReactive(dgeObj(), {
    calcNormFactors(dgeObj())
  })

  ## compute the log-CPM of the counts for unnormalized and normalized objects
  logCpm <- eventReactive(dgeObj(), {
    cpm(dgeObj(), log = TRUE)
  })

  logCpmNorm <- eventReactive(dgeObjNorm(), {
    cpm(dgeObjNorm(), log = TRUE)
  })

  # Plot size and logCpm distribution
  ## size of libraries
  output$libSizes <- renderPlot({
    barplot(
      dgeObj()$samples$lib.size,
      names = colnames(dgeObj()),
      las = 2,
      main = "Barplot of library sizes"
    )
  })

  ## log-CPM distribution of libraries
  output$logCpmNormDistrib <- renderPlot({
    boxplot(
      logCpmNorm(),
      xlab = "",
      ylab = "Log2 counts per million",
      las = 2,
      main = "Boxplots of logCpms (normalised)"
    )
    abline(h = median(logCpmNorm()), col = "blue")
  })

  # Return results
  results <- list(dgeObj, dgeObjNorm, logCpm, logCpmNorm)
  names(results) <- c("dgeObj", "dgeObjNorm", "logCpm", "logCpmNorm")
  return(results)
}