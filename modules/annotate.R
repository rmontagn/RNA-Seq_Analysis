### -------------------------------------------------------------------- ###
###                         MODULE ANNOTATE
### -------------------------------------------------------------------- ###
library(edgeR)

### FUNCTIONS
### -------------------------------------------------------------------- ###
annotateGenes <- function(df) {
  # Set up connection to ensembl database
  usedMart=useMart("ENSEMBL_MART_ENSEMBL")
  # list the available datasets (species)
  # listDatasets(ensembl) %>%
  #   filter(str_detect(description, "Human"))
  
  # Specify a data set to use
  ensembl = useDataset("hsapiens_gene_ensembl", mart=usedMart)
  
  # Set the filter type and values
  filterType <- "ensembl_gene_id"
  filterValues <- gsub("\\..+", "", rownames(df))
  
  # Check the available "attributes" - things you can retreive
  # listAttributes(ensembl)[,c(1,2)] %>%
  # head(20)
  
  # Set the list of attributes
  attributeNames <- c('ensembl_gene_id', 'external_gene_name')#'ensembl_transcript_id', 'go_id', 'external_gene_name', 'description')
  
  # Run the query
  annot <- getBM(attributes=attributeNames,
                 filters = filterType,
                 values = filterValues,
                 mart = ensembl)
  
  return(annot)
}

### UI
### -------------------------------------------------------------------- ###
genesAnnotationInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  return(NULL)
}


### SERVER
### -------------------------------------------------------------------- ###
genesAnnotation <- function(input, output, session, button, condition, counts) {
  
  annotatedGenes <- eventReactive(button(), {
    if(condition()) {
    print("Sending Query to Biomart")
    annotateGenes(counts())
    } else {
      NULL
    }
  })
  
  return(annotatedGenes)
}