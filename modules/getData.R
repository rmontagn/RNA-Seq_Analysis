### -------------------------------------------------------------------- ###
###                     MODULE GET DATA
### -------------------------------------------------------------------- ###
library(edgeR)

### FUNCTIONS
### -------------------------------------------------------------------- ###
# Get raw counts table
getRawData <- function(countsFile) {
  countsData <-
    read.delim(
      as.character(countsFile),
      stringsAsFactors = FALSE,
      header = T,
      sep = "\t",
    )
  rownames(countsData) <- countsData[, 1]
  countsData <- countsData[, -1]
  return(as.data.frame(countsData))
}

# Remove duplicates from the counts table
removeDuplicates <- function(countsData) {
  # genes <- gsub("\\..+", "", rownames(countsData))
  # isDup <- duplicated(genes)
  # dup <- results$genes[isDup]
  # results[results$genes%in%dup,]
  countsData$gene <- gsub("\\..+", "", rownames(countsData))
  countsData <- unique(countsData)
  toRemove <- grep('PAR_Y', rownames(countsData))
  if(!is_empty(toRemove)) {
    print(c("to remove: ", toRemove))
    countsData <- countsData[-toRemove,]
  }
  rownames(countsData) <- countsData$gene
  countsData <- countsData[,-ncol(countsData)]
  return(countsData)
}

# Remove lowly expressed genes
removeLowlyExpressedGenes <- function(countsData) {
  ## here genes that have 4 counts or less in at least 5% of patients
  nb.patients <- as.integer(dim(countsData)[2])
  thresh <- countsData > 4
  ## rowSums(thresh)
  keep <- rowSums(thresh) >= floor(0.05 * nb.patients)
  counts.keep <- countsData[keep, ]
  
  return(counts.keep)
}

makeFeaturesList <- function(df) {
  chlst <- c("", colnames(df))
  names(chlst) <- c("", colnames(df))
  return(chlst)
}


### UI
### -------------------------------------------------------------------- ###
getModuleDataInput <- function(id) {
  
  ns <- NS(id)
  
  tagList(
      # Input for the counts table
      fileInput(ns("countsFile"),
                label = "Select a RNA-seq count table \n "),
      # tags$param("Samples must be columns and genes must be rows. The first line must be a header with the sample names. The first columns must contain the gene identifiers."),
      
      # Input for the features table
      fileInput(ns("featuresFile"),
                label = "Select the associated clinical features"),
      
      # Checkbox for early annotation
      checkboxInput(ns("annotateGenes"), 
                    label = "annotated genes (may take some time)",
                    value = FALSE),
      
      # Action button
      actionButton(ns("explore"), label = "Exploratory Analysis")
  )
}


### SERVER
### -------------------------------------------------------------------- ###
getModuleData <- function(input, output, session) {
  
  # The properties of the input count table
  # file <- reactive(input$featuresFile)
  observeEvent(input$explore, {
    print(input$countsFile$datapath)
  })

  # The counts table
  counts <- eventReactive(input$explore, {
    print("get the counts")
    data.matrix(removeLowlyExpressedGenes(removeDuplicates(getRawData(input$countsFile$datapath))))
  })  
 
  # The features table
  features <- eventReactive(input$explore, {
    print("get the features")
    as.data.frame(getRawData(input$featuresFile$datapath))
  })
  
  # The list of clinical features
  chlst <- eventReactive(input$explore, {
    print("make list of features")
    makeFeaturesList(as.data.frame(features()))
  })
  
  results <- list(counts, features, chlst, reactive(input$explore), reactive(input$annotateGenes))
  names(results) <- c("counts", "features", "chlst", "inputExplore", "inputAnnotatedGenes")
  return(results)
  
}