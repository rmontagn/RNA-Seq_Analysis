### -------------------------------------------------------------------- ###
###                            MODULE HEATMAP
### -------------------------------------------------------------------- ###
# Compute the most variable gene in the dge object

library(edgeR)
library(RColorBrewer)
rdybu.palette <- brewer.pal(11, "RdYlBu")
heatmap.colors <- colorRampPalette(rdybu.palette)

### FUNCTIONS
### -------------------------------------------------------------------- ###


### UI
### -------------------------------------------------------------------- ###
mostVariableGenesOutput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)

  plotOutput(ns("heatmap"))

}


### SERVER
### -------------------------------------------------------------------- ###
mostVariableGenes <- function(input, output, session, logCpmNorm) {

  # Output the Heatmap of the 500 most variable genes
  output$heatmap <- renderPlot({
    ## get the 500 most variable genes
    genesVariance <- apply(logCpmNorm(), 1, var)
    mostVariableGenes <- names(
      sort(
        genesVariance, decreasing = TRUE
        )
      )[1:500]

    ## get the logCpmNorm of these genes
    mostVariableLcpm <- logCpmNorm()[mostVariableGenes, ]

    ## plot their heatmap
    heatmap.2(
      mostVariableLcpm, 
      col = rev(heatmap.colors(50)),
      trace = "none",
      main = "Top 500 most variable genes across samples",
      # ColSideColors=col.cell,
      scale = "row"
      )
    })
}