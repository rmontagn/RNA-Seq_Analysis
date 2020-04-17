### -------------------------------------------------------------------- ###
###                    MODULE FUNCTION ENRICHMENT 
### -------------------------------------------------------------------- ###


### FUNCTIONS
### -------------------------------------------------------------------- ###
# Find the enriched functions for the de genes
getEnrichedFunctions <- function(df) {
  genes <- df$FDR < 0.01
  names(genes) <- row.names(df)
  pwf <- nullp(genes, "hg19", "ensGene")
  go.results <- goseq(pwf, "hg19","ensGene", test.cats = c("GO:BP"))
  enriched.GO <- go.results$category[p.adjust(go.results$over_represented_pvalue,method="BH")<.01]
  enrichedGoResults <- as.data.frame(merge(go.results, as.data.frame(enriched.GO), by.x="category", by.y="enriched.GO"))

  return (enrichedGoResults)
}

### UI
### -------------------------------------------------------------------- ###
functionEnrichmentOutput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  fluidRow(plotOutput(ns("deEnrichedFunctions")))
}


### SERVER
### -------------------------------------------------------------------- ###
functionEnrichment <- function(input, output, session, button, dgeDf) {
  
  # Output the function enrichment
  observeEvent(button(), {
    output$deEnrichedFunctions <- renderPlot({
      getEnrichedFunctions(dgeDf) %>%
        top_n(10, wt=-over_represented_pvalue) %>%
        mutate(hitsPerc=numDEInCat*100/numInCat) %>%
          ggplot(aes(x=hitsPerc, 
                     y=term, 
                     colour=over_represented_pvalue, 
                     size=numDEInCat)) +
          geom_point() +
          expand_limits(x=0) +
          labs(x="Hits (%)", y="GO term", colour="p value", size="Count")
    })
  })
}
