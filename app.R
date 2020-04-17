# ==============================
# Load packages
# ==============================
list.of.packages <- c("ggplot2", "Rcpp", "rmarkdown", "gplots", "RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)

list.of.bioconductor.packages <-
  c("limma",
    "Glimma",
    "edgeR",
    "DESeq",
    "org.Hs.eg.db",
    "Homo.sapiens",
    "GO.db")
new.packages <- list.of.bioconductor.packages[!(list.of.bioconductor.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  BiocManager::install(new.packages)

rm(list = ls())

library(edgeR)
library(Glimma)
library(biomaRt)
library(goseq)
# library(GO.db)

library(tidyverse)
library(gplots)
library(ggplot2)
library(RColorBrewer)

library(shiny)
source("modules/getData.R")
source("modules/annotate.R")
source("modules/dgeObj.R")
source("modules/mds.R")
source("modules/heatmap.R")
source("modules/diffExp.R")
source("modules/smear.R")
source("modules/volcano.R")
source("modules/enrichment.R")


### -------------------------------------------------------------------- ###
###                            UI
### -------------------------------------------------------------------- ###
ui <- fluidPage(
  
  titlePanel("RNA-seq data Analyzer"),

  sidebarLayout(
    ### Sidebar
    ### -------------------------------------------------------------------- ###
    sidebarPanel(
      # Subpanel1: Input the data
      # --------------------------------------
      tags$h4(style="text-align: center", "Load data"),
      inputPanel(getModuleDataInput("inputData")),
      genesAnnotationInput("annot"), # just here to make module annotate.R work
      br(),
    
      # Subpanel2: Choose the feature and condition to test against the other group(s)
      # --------------------------------------
      tags$h4(style="text-align: center", "Differential Expression"),
      wellPanel(
        # Clinical feature for patient groups
        selectInput(inputId = "mdsGroupingFeature", label = "Choose a clinical feature to group your samples", choices = ""),
        # Radio button to choose condition to test
        conditionalPanel(condition = "input.mdsGroupingFeature != ''",radioButtons(inputId = "DeCondition", label = "Choose a condition to test against the others", choices = "")),
      
        # Action buttons
        actionButton(inputId = "runDe", label = "Differential Expression Analysis"),
        actionButton(inputId = "save.de", label = "Save model"),

        # Input for loading data saved as a rds file
        fileInput(inputId = "loadDE", label = "Load model (rds file)")
      ),
      
      # Subpanel3: trigger function enrichment analysis
      # --------------------------------------
      tags$h4(style="text-align: center", "Functional Enrichment"),
      wellPanel(
        conditionalPanel(condition = "input.runDe > 0", actionButton("enrichment", label="Function Enrichment"))
      )
    ), 
  
    ### Main Panel
    ### -------------------------------------------------------------------- ###
    mainPanel(
      # Subpanel1: Explore the data
      # --------------------------------------
      wellPanel(
        "Data Exploration",
        # First row: simple size and logCpm distributions
        computeDgeObjOutput("dgeObjs"),
        br(),
        # Second row: MDS and most variable genes heatmap
        fluidRow(
          column(6, align = "center", mdsOutput("mds")),
          column(6, mostVariableGenesOutput("heatmap")),
        )
      ),
      
      # Subpanel2: Differential expression
      # --------------------------------------
      wellPanel(
        "Differential Gene Expression",
        conditionalPanel(
          condition = "input.runDe > 0",
          fluidRow(computeDiffExpOutput("diffExp")),
          br(),
          fluidRow(
            column(6, computePlotSmearOutput("smear")),
            column(6, computeVolcanoPlotOutput("volcano"))
          ),
        )
      ),
      
      # Subpanel3: Functional enrichment
      # --------------------------------------
      wellPanel(
        "Functional Enrichment",
        conditionalPanel(
          condition = "input.enrichment > 0",
          functionEnrichmentOutput("enrichedFunctions")
        )
      )
    )
  )
)
                

### -------------------------------------------------------------------- ###
###                            SERVER
### -------------------------------------------------------------------- ###
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 1000 * 1024 ^ 2)
  ## Get some nice colours using Rcolorbrewer
  rdybu.palette <- brewer.pal(11, "RdYlBu")
  set2.palette <- brewer.pal(8, "Set2")
  mds.colors <- colorRampPalette(set2.palette)
  heatmap.colors <- colorRampPalette(rdybu.palette)
  
  ### Get the data
  ### -------------------------------------------------------------------- ###
  results <- callModule(getModuleData, "inputData")
  counts <- results$counts
  features <- results$features
  chlst <- results$chlst
  inputExplore <- results$inputExplore
  inputAnnotatedGenes <- results$inputAnnotatedGenes
  
  ### Annotate the gene id
  ### -------------------------------------------------------------------- ###
  annotatedGenes <- callModule(genesAnnotation, "annot", inputExplore, inputAnnotatedGenes, counts)

  ### Get DGE Object
  ### -------------------------------------------------------------------- ###
  dgeObjs <- callModule(computeDgeObj, "dgeObjs", inputExplore, counts(), features(), annotatedGenes())
  dgeObj <- dgeObjs$dgeObj
  dgeObjNorm <- dgeObjs$dgeObjNorm
  logCpm <- dgeObjs$logCpm
  logCpmNorm <- dgeObjs$logCpmNorm
  
  ### Compute MDS
  ### -------------------------------------------------------------------- ###
  # Update the list of possible grouping features for MDS (to refactor ?)
  observeEvent(chlst(), {
    updateSelectInput(session, "mdsGroupingFeature", choices = chlst())
  })
  
  # Get the grouping feature for MDS and next differential expression analysis
  mdsGroupingFeature <- reactive(input$mdsGroupingFeature)
  
  # Limma plain MDS
  callModule(mds, "mds", mdsGroupingFeature=reactive(input$mdsGroupingFeature), dgeObjNorm=dgeObjNorm)
  
  
  ### Compute 500 most variable genes
  ### -------------------------------------------------------------------- ###
  callModule(mostVariableGenes, "heatmap", logCpmNorm)
  
  
  ### Compute differential expression
  ### -------------------------------------------------------------------- ###
  # Update the list of possible conditions for the selected grouping feature (to refactor ?)
  observeEvent(mdsGroupingFeature(), {
    updateRadioButtons(session, "DeCondition", choices = unique(dgeObj()$samples[, as.character(input$mdsGroupingFeature)]))
  })
  
  # Compute differential expression, plot the summary of results and the histograms of p-values
  rv <- callModule(computeDiffExp, "diffExp", reactive(input$runDe), dgeObjNorm, mdsGroupingFeature, reactive(input$DeCondition))
  
  # Save the lrt model
  observeEvent(input$save.de, {
    lrt <- lrt()
    saveRDS(lrt, file="./model.rds")
  })
  
  
  ### Compute smear and volcano plot
  ### -------------------------------------------------------------------- ###
  callModule(computePlotSmear, "smear", reactive(rv$lrtModel), reactive(rv$dgeChangedGenes), counts, reactive(rv$dgeDecide), dgeObjNorm, reactive(input$mdsGroupingFeature))
  callModule(computeVolcanoPlot, "volcano", reactive(rv$dgeDf), reactive(rv$dgeChangedGenes))
  
  ### Compute functions enrichment
  ### -------------------------------------------------------------------- ###
  callModule(functionEnrichment, "enrichedFunctions", reactive(input$enrichment), rv$dgeDf)
}

shinyApp(ui = ui, server = server)