### -------------------------------------------------------------------- ###
###                            MODULE DGEOBJ
### -------------------------------------------------------------------- ###

library(edgeR)

### FUNCTIONS
### -------------------------------------------------------------------- ###


### UI
### -------------------------------------------------------------------- ###
computeDgeObjOutput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
return(NULL)
}


### SERVER
### -------------------------------------------------------------------- ###
computeDgeObj <- function(input, output, session, button, counts, features, annot) {
  
  # Create dge object (an EdgeR object) from a RNA-seq matrix
  dgeObj <- eventReactive(button(), {
    if(is.null(annot)){
      DGEList(counts[-1,], samples = features)
    } else {
      DGEList(counts[-1,], samples = features, genes = annot[,2])
    }
  })
  
  # Normalize counts  
  dgeObjNorm <- eventReactive(dgeObj(), {
   calcNormFactors(dgeObj())
  })
  
  # Compute the log-CPM of the counts for each object
  logCpm <- eventReactive(dgeObj(), {
    cpm(dgeObj(), log = TRUE)
  })

  logCpmNorm <- eventReactive(dgeObjNorm(), {
    cpm(dgeObjNorm(), log = TRUE)
  })
  
  # Return results
  results <- list(dgeObj, dgeObjNorm, logCpm, logCpmNorm)
  names(results) <- c("dgeObj", "dgeObjNorm", "logCpm", "logCpmNorm")
  return(results)
}