### -------------------------------------------------------------------- ###
###                            MODULE SMEAR
### -------------------------------------------------------------------- ###
library(edgeR)

### FUNCTIONS
### -------------------------------------------------------------------- ###


### UI
### -------------------------------------------------------------------- ###
computeVolcanoPlotOutput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  tagList(
    fluidRow(plotOutput(ns("deVolcano"))),
    br(),
    fluidRow(actionButton(ns("runGlimmaVolcano"), label = "Glimma Volcano Plot")),
  )
}


### SERVER
### -------------------------------------------------------------------- ###
computeVolcanoPlot <- function(input, output, session, dgeDf, dgeChangedGenes) {

  output$deVolcano <- renderPlot({
    if(!is.null(dgeDf())){
      signif <- -log10(dgeDf()$FDR)
      plot(
        dgeDf()$logFC,
        signif,
        pch = 16,
        main = "Volcano plot",
        xlab = "logFC",
        ylab = expression("-log"["10"] * "(FDR)")
      )
      points(dgeDf()[dgeChangedGenes(), "logFC"], -log10(dgeDf()[dgeChangedGenes(), "FDR"]), pch = 16, col = "red")
      abline(v = -2, col = "blue")
      abline(v = 2, col = "blue")
      abline(h = 1.3, col = "green")
    }
  })
}