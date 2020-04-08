### -------------------------------------------------------------------- ###
###                            MODULE VOLCANO
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
    htmlOutput(ns("glimmaVolcano"))
  )
}


### SERVER
### -------------------------------------------------------------------- ###
computeVolcanoPlot <- function(input, output, session, dgeDf, dgeChangedGenes) {

  output$deVolcano <- renderPlot({
    if(!is.null(dgeDf())){
      colors <- as.numeric(dgeDf()$FDR <= 0.05 & dgeDf()$logFC > 0) + 2*as.numeric(dgeDf()$FDR <= 0.05 & dgeDf()$logFC < 0) + 1
      signif <- -log10(dgeDf()$FDR)
      plot(
        dgeDf()$logFC,
        signif,
        pch = 16,
        main = "Volcano plot",
        xlab = "logFC",
        ylab = expression("-log"["10"] * "(FDR)"),
        col = c("black", "red", "blue")[colors]
      )
      abline(v = -2, col = "blue")
      abline(v = 2, col = "blue")
      abline(h = 1.3, col = "green")
    }
  })
  
  observeEvent(input$runGlimmaVolcano, {
    output$glimmaVolcano <- renderUI({
      dgeStatus <- as.numeric(dgeDf()$FDR <= 0.05 & dgeDf()$logFC > 0) - as.numeric(dgeDf()$FDR <= 0.05 & dgeDf()$logFC < 0)
      ga2=data.frame(GeneID=rownames(dgeDf()), rownames=rownames(dgeDf()))
      glXYPlot(x= dgeDf()$logFC,
               y=-log10(dgeDf()$FDR),
               xlab="logFC", 
               ylab="-log(FDR)", 
               status=as.numeric(dgeStatus),
               anno=ga2)
    })
  })
}