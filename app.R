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

# ============================================================
# SERVER FUNCTIONS
# ============================================================

# # Get raw counts table
# getRawData <- function(countsFile) {
#   countsData <-
#     read.delim(
#       as.character(countsFile),
#       stringsAsFactors = FALSE,
#       header = T,
#       sep = "\t",
#     )
#   # colnames(seqdata) <- seqdata[1, ]
#   rownames(countsData) <- counts.data[, 1]
#   countsData <- countsData[, -1]
#   # counts.data <- data.matrix(counts.data)
#   return(as.data.frame(countsData))
# }

# annotateGenes <- function(df) {
#   # Set up connection to ensembl database
#   usedMart=useMart("ENSEMBL_MART_ENSEMBL")
#   # list the available datasets (species)
#   # listDatasets(ensembl) %>%
#   #   filter(str_detect(description, "Human"))
# 
#   # Specify a data set to use
#   ensembl = useDataset("hsapiens_gene_ensembl", mart=usedMart)
# 
#   # Set the filter type and values
#   filterType <- "ensembl_gene_id"
#   filterValues <- gsub("\\..+", "", rownames(df))
# 
#   # Check the available "attributes" - things you can retreive
#   # listAttributes(ensembl)[,c(1,2)] %>%
#   #   head(20)
# 
#   # Set the list of attributes
#   attributeNames <- c('ensembl_gene_id', 'external_gene_name')#'ensembl_transcript_id', 'go_id', 'external_gene_name', 'description')
# 
#   # Run the query
#   annot <- getBM(attributes=attributeNames,
#                  filters = filterType,
#                  values = filterValues,
#                  mart = ensembl)
#   
#   return(annot)
# }

# removeDuplicates <- function(counts.data) {
#   # genes <- gsub("\\..+", "", rownames(counts.data))
#   # isDup <- duplicated(genes)
#   # dup <- results$genes[isDup]
#   # results[results$genes%in%dup,]
#   counts.data$gene <- gsub("\\..+", "", rownames(counts.data))
#   counts.data <- unique(counts.data)
#   toRemove <- grep('PAR_Y', rownames(counts.data))
#   if(!is_empty(toRemove)) {
#     print(c("to remove: ", toRemove))
#     counts.data <- counts.data[-toRemove,]
#   }
#   rownames(counts.data) <- counts.data$gene
#   counts.data <- counts.data[,-ncol(counts.data)]
#   return(counts.data)
# }

# Remove lowly expressed genes
# remove.lowly.expressed.genes <- function(counts.data) {
#   ## remove lowly expressed genes
#   ## here genes that have 4 counts or less in at least 5% of patients
#   nb.patients <- as.integer(dim(counts.data)[2])
#   thresh <- counts.data > 4
#   # rowSums(thresh)
#   keep <- rowSums(thresh) >= floor(0.05 * nb.patients)
#   # summary(keep)
#   counts.keep <- counts.data[keep, ]
#   # logcounts.keep <- log(counts.keep)
#   
#   return(counts.keep)
# }

# # Compute log-CPM from a dge object
# log.cpm.distrib <- function(dgeObj) {
#   return(cpm(dgeObj, log = TRUE))
# }

# Turns a vector of character into a list that can be displayed in a selectInput object
  # make.features.list <- function(df) {
  #   chlst <- c("", colnames(df))
  #   names(chlst) <- c("", colnames(df))
  #   return(chlst)
  # }

# Deduce LRT model from dge object
# analyze.de <- function(dgeObj, feature, condition) {
#   dgeObj <- estimateCommonDisp(dgeObj)
#   dgeObj <- estimateGLMTrendedDisp(dgeObj)
#   dgeObj <- estimateTagwiseDisp(dgeObj)
#   # plotBCV(dgeObj)
#   feature.vector <- dgeObj$samples[, as.character(feature)]
#   feature.vector <- feature.vector == condition
#   design <- model.matrix(~ feature.vector)
#   fit <- glmFit(dgeObj, design)
#   lrt <- glmLRT(fit, coef = 2)
#   return(lrt)
# }

# Find the enriched functions for the de genes
getEnrichedFunctions <- function(df) {
  observeEvent(df(),{
    df <- isolate(df())
    genes <- df$FDR < 0.01
    names(genes) <- row.names(df)
    pwf <- nullp(genes, "hg19", "ensGene")
    go.results <- goseq(pwf, "hg19","ensGene", test.cats = c("GO:BP"))
    print(length(p.adjust(go.results$over_represented_pvalue,method="BH")))
    enriched.GO <- go.results$category[p.adjust(go.results$over_represented_pvalue,method="BH")<.01]
    enrichedGoResults <- as.data.frame(merge(go.results, as.data.frame(enriched.GO), by.x="category", by.y="enriched.GO"))
    getEnrichedFunctions() %>%
      top_n(10, wt=-over_represented_pvalue) %>%
      mutate(hitsPerc=numDEInCat*100/numInCat) %>% test
  
    return(test)
  })
}

# Plot the enriched functions with ggplot2
# plotEnrichedFunctions <- function(enrichedGoResults) {
#   plot <- enrichedGoResults %>%
#     top_n(10, wt=-over_represented_pvalue) %>%
#     mutate(hitsPerc=numDEInCat*100/numInCat) %>%
#     ggplot(aes(x=hitsPerc,
#                y=term,
#                colour=over_represented_pvalue,
#                size=numDEInCat)) +
#     geom_point() +
#     expand_limits(x=0) +
#     labs(x="Hits (%)", y="GO term", colour="p value", size="Count")
#   # ggsave(file.path(directory, "most_enriched_functions.png"))
#   return(plot)
# }


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
          column(6, align = "center", fluidRow(mdsOutput("mds"))),
                 # fluidRow(
                 #   conditionalPanel(
                 #     condition = "input.explore > 0",
                 #     actionButton(inputId = "run.glimma", label = "Glimma MDS")
                 #   )
          column(6, mostVariableGenesOutput("heatmap")),
        # htmlOutput("glimma") # , inline=TRUE)
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
        #   column(6,
        #        fluidRow(plotOutput(outputId = "de.volcano")),
        #        fluidRow(actionButton(inputId = "runGlimmaVolcano", label = "Glimma Volcano Plot"))
        #   )
        # ),
        # br(),
        # fluidRow(
        #   plotOutput(outputId = "deEnrichedFunctions")
        # ),
        # htmlOutput(outputId="de.volcano")),
        # htmlOutput("glimmaSmear"),
        # htmlOutput("glimmaVolcano"),
        # tableOutput(outputId = "enrichedFunctions")#,
        fluidRow(plotOutput(outputId = "deEnrichedFunctions"))
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
  
  # The properties of the input count table
  #   file <- reactive(input$features.file)
  #   
  #   # The counts table
  #   counts <- eventReactive(input$explore, {
  #     data.matrix(remove.lowly.expressed.genes(removeDuplicates(getRawData(input$countsFile$datapath))))
  #   })
  # 
  # # The features table
  # features <- eventReactive(input$explore, {
  #   as.data.frame(getRawData(input$features.file$datapath))
  # })
  
  # # The list of clinical features
  # chlst <- eventReactive(input$explore, {
  #   make.features.list(as.data.frame(features()))
  # })
  
  # annotatedGenes <- eventReactive(inputExplore(), {
  #   if(inputAnnotatedGenes()) {
  #     print("Sending Query to Biomart")
  #     annotateGenes(counts())
  #   }
  #   else {
  #     NULL
  #   }
  # })
  
  # keep <- eventReactive(input$explore, {
  #   intersect(names(counts()), names(t(features())))
  # })
  
  # The dge obect created by EdgeR after loading the data
  # and removing lowly expressed genes

  # dgeObj <- eventReactive(input$explore, {
  #   create.dge.object(as.matrix(counts()), features(), annotatedGenes())
  # })
  # 
  # dgeObjNorm <- eventReactive(dgeObj(), {
  #   calcNormFactors(dgeObj())
  # })
  
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
  
  # # Glimma interactive MDS
  # observeEvent(input$run.glimma, {
  #   if (input$mdsGroupingFeature != "") {
  #     glMDSPlot(
  #       dgeObjNorm(),
  #       labels = rownames(dgeObjNorm()$samples),
  #       groups = dgeObjNorm()$samples[, as.character(input$mdsGroupingFeature)]
  #     )
  #   } else {
  #     glMDSPlot(dgeObjNorm(), labels = rownames(dgeObjNorm()$samples))
  #   }
  # })
  # 
  # output$glimma <- renderUI({
  #   includeHTML("glimma-plots/MDS-Plot.html")
  # })
  
  
  ### Compute 500 most variable genes
  ### -------------------------------------------------------------------- ###
  callModule(mostVariableGenes, "heatmap", logCpmNorm)
  
  # logCpm <- eventReactive(dgeObj(), {
  #   cpm(dgeObj(), log = TRUE)
  # })
  # 
  # logCpmNorm <- eventReactive(dgeObjNorm(), {
  #   cpm(dgeObjNorm(), log = TRUE)
  # })
  # limma.groups <- reactive({
  #   groups <- as.factor(merge(dge()$sample, features(),by="row.names")[,as.character(input$mdsGroupingFeature)])
  # })
  #

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
  
  # lrtModel <- dgeResults$lrtModel
  # dgeDf <- dgeResults$dgeDf
  # dgeDecide <- dgeResults$dgeDecide
  # lrt <- eventReactive(input$run.de, {
  #   analyze.de(dgeObjNorm(),
  #              input$mdsGroupingFeature,
  #              as.character(input$DE.condition))
  # })
  # 
  # de.df <- eventReactive(lrt(), {
  #   as.data.frame(topTags(lrt(), n = Inf))
  # })
  # 
  # deDecideTest <- eventReactive(de.df(), {
  #   decideTestsDGE(lrt(), p = 0.05, adjust = "BH")
  # })
  # 
  # deChangedGenes <-
  #   eventReactive(deDecideTest(), {
  #     rownames(dgeObjNorm()[as.logical(deDecideTest()), ])
  #   })

  # getAnnotatedGenes <- function(my_keys, unannotatedDf) {
  #   rownames(unannotatedDf) <- gsub("\\..*","", rownames(unannotatedDf))
  #   # Problem: there are duplicate rownames
  #   ann <- as.data.frame(select(org.Hs.eg.db, keys=my_keys, keytype="ENSEMBL", columns=c("ENSEMBL","SYMBOL","GENENAME", "GO")))
  #   print(head(ann))
  #   ann <- ann[ann$ONTOLOGY == "BP",]
  #   library(dplyr)
  #   ann <- distinct(ann, ENSEMBL, .keep_all= TRUE)
  #   detach("package:dplyr", unload=TRUE)
  #   results.annotated <- merge(ann, unannotatedDf, by.x="ENSEMBL", by.y="ENSEMBL", all.x=T, all.y=T)
  #   results.annotated <- results.annotated[-nrow(results.annotated),]
  #   print(head(results.annotated))
  #   return(results.annotated)
  # }


  # observeEvent(input$save.de, {
  #   lrt <- lrt()
  #   saveRDS(lrt, file="./model.rds")
  # })
  
  # lrt <- eventReactive(input$loadDE, {
  #   readRDS(lrt, file=input$loadDE$datapath)
  # })
  
    # enrichedFunctions <- eventReactive(deDecideTest(), {
    #   getEnrichedFunctions(deDecideTest())
    # })

  # PLOTS
  # ============================================================
  
  # Panel Explore
  # --------------------------------------
  # Plot the counts table file propertiess
  # output$file <- renderPrint(file())

  ## MDS
  ## compute and plot the MDS
  # output$selected.mds.group <- renderPrint({
  #   as.numeric(as.factor(features()[,as.character       (input$mdsGroupingFeature)]))
  # })
  # output$mds <- renderPlot({
  #   if (input$mdsGroupingFeature != "") {
  #     # plotMDS(dge(), labels=labels, groups=group, folder="mds")
  #     feat <-
  #       as.factor(dgeObjNorm()$samples[, as.character(input$mdsGroupingFeature)])
  #     plotMDS(dgeObjNorm(),
  #             pch = 16,
  #             main = "MDS",
  #             col = set2.palette[as.numeric(feat)])
  #     legend(
  #       "topleft",
  #       col = set2.palette,
  #       pch = 16,
  #       legend = levels(as.factor(feat))
  #     )
  #   } else {
  #     plotMDS(dgeObjNorm(), pch = 16, main = "MDS")
  #   }
  # })
  
  #-- Panel Differential Expression ---------------------------------------
  
  ## update the possible conditions to choose for DE
  # observeEvent(mdsGroupingFeature(), {
  #   updateRadioButtons(session, "DE.control", choices=unique(dgeObj()$samples[,as.character(input$mdsGroupingFeature)]))
  # })
  
  # ## display the summary of the DE analysis
  # output$de.summary <- renderPrint({
  #   summary(deDecideTest())
  # })
  
  # Plot the histogram of p-values for DE genes
  # output$de.pval <- renderPlot({
  #   # print(head(de.df()))
  #   hist(de.df()$PValue,
  #        las = 2,
  #        main = "Histogram of p-values",
  #        xlab = "p-values")
  # })
  # 
  # ## plot the histogram of FDR
  # output$de.fdr <- renderPlot({
  #   hist(de.df()$FDR,
  #        las = 2,
  #        main = "Histogram of FDR",
  #        xlab = "FDR")
  #   abline(v = 0.05, col = "blue")
  #   
  # })
  
  # Plot the smear with selected genes
  # output$de.smear <- renderPlot({
  #   plotSmear(dgeResults$lrtModel, de.tags = dgeResults$dgeDecide, main = "Smear plot")
  #   # points(de.changed.genes()$logCpm, de.changed.genes()$logFC, col="red")#, labels = significant$SYMBOL,col="red")
  # })
  # observeEvent(input$runGlimmaSmear, {
  #   output$glimmaSmear <- renderUI({
  #     glMDPlot(lrt(),
  #              counts=counts()[-1,],
  #              xlab="logFC", 
  #              ylab="-log(FDR)", 
  #              status=deDecideTest(), 
  #              anno=annotatedGenes(), 
  #              groups=as.factor(dgeObjNorm()$samples[, as.character(input$mdsGroupingFeature)]),
  #              side.main="external_gene_nama"
  #              )
  #   })
  # })

  # Plot the volcano plot with selected genes
  # output$de.volcano <- renderPlot({
  #   signif <- -log10(de.df()$FDR)
  #   plot(
  #     de.df()$logFC,
  #     signif,
  #     pch = 16,
  #     main = "Volcano plot",
  #     xlab = "logFC",
  #     ylab = expression("-log"["10"] * "(FDR)")
  #   )
  #   points(de.df()[deChangedGenes(), "logFC"], -log10(de.df()[deChangedGenes(), "FDR"]), pch = 16, col = "red")
  #   abline(v = -2, col = "blue")
  #   abline(v = 2, col = "blue")
  #   abline(h = 1.3, col = "green")
  # })

  # observeEvent(input$runGlimmaVolcano, {
  #   output$glimmaVolcano <- renderUI({  
  #     ga2=data.frame(GeneID=rownames(de.df()), rownames=rownames(de.df()))
  #     glXYPlot(x= de.df()$logFC,
  #              y=-log10(de.df()$FDR),
  #              xlab="logFC", 
  #              ylab="-log(FDR)", 
  #              status=as.numeric(de.df()$FDR <= 0.05), 
  #              anno=ga2)
  #   })
    
    # glXYPlot(df$logFC,-log10(df$FDR),
    #          xlab="logFC",
    #          ylab="-log(FDR)",
    #          status=as.numeric(df$FDR <= 0.05),
    #          anno=ga2)
  # })
  # output$annotatedGenes <- renderTable(head(deDf()))

  output$deEnrichedFunctions <- renderPlot({
    list(ggplot(getEnrichedFunctions(reactive(rv$dgeDf)), aes_string(x=hitsPerc,
                 y=term,
                 colour=over_represented_pvalue,
                 size=numDEInCat)
           +geom_point()
           + expand_limits(x=0)
           + labs(x="Hits (%)", y="GO term", colour="p value", size="Count")
    ))
  })
}

shinyApp(ui = ui, server = server)

# detags <- rownames(dgeObj)[as.logical(dt)]
# ann <- as.data.frame(select(org.Mm.eg.db, keys=detags, keytype="ENSEMBL", columns=c("ENSEMBL","SYMBOL","GENENAME", "GO")))
# ann <- ann[ann$ONTOLOGY == "BP",]
# library(dplyr)
# ann <- distinct(ann, ENSEMBL, .keep_all= TRUE)
# detach("package:dplyr", unload=TRUE)
# results.annotated <- merge(ann, results, by.x="ENSEMBL", by.y="ENSEMBL", all.x=T, all.y=T)
# results.annotated <- results.annotated[-nrow(results.annotated),]