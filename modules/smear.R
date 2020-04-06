### -------------------------------------------------------------------- ###
###                            MODULE SMEAR
### -------------------------------------------------------------------- ###

library(edgeR)

### FUNCTIONS
### -------------------------------------------------------------------- ###


### UI
### -------------------------------------------------------------------- ###
computePlotSmearOutput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  tagList(
    fluidRow(plotOutput(ns("deSmear"))),
    br(),
    fluidRow(actionButton(ns("runGlimmaSmear"), label = "Glimma Smear")),
    htmlOutput(ns("glimmaSmear"))
  )
}
  
  
### SERVER
### -------------------------------------------------------------------- ###
computePlotSmear <- function(input, output, session, lrtModel, dgeChangedGenes, counts, dgeDecide, dgeObjNorm, mdsGroupingFeature) {

  output$deSmear <- renderPlot({
    # print(dgeChangedGenes())
    if(!is.null(lrtModel())){
    plotSmear(lrtModel(), main = "Smear plot", de.tags = dgeChangedGenes())
    }
  })
  
  observeEvent(input$runGlimmaSmear, {
    output$glimmaSmear <- renderUI({
      glMDPlot(lrtModel(),
               counts=counts()[-1,],
               xlab="logFC",
               ylab="-log(FDR)",
               status=dgeDecide(),
               # anno=annotatedGenes(),
               groups=as.factor(dgeObjNorm()$samples[, as.character(mdsGroupingFeature())]),
               side.main="external_gene_nama"
      )
    })
  })

}