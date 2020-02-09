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


# ============================================================
# SERVER FUNCTIONS
# ============================================================

# Get raw counts table
getRawData <- function(counts.file) {
  counts.data <-
    read.delim(
      as.character(counts.file),
      stringsAsFactors = FALSE,
      header = T,
      sep = "\t",
    )
  # colnames(seqdata) <- seqdata[1, ]
  rownames(counts.data) <- counts.data[, 1]
  counts.data <- counts.data[, -1]
  # counts.data <- data.matrix(counts.data)
  return(as.data.frame(counts.data))
}

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
  #   head(20)

  # Set the list of attributes
  attributeNames <- c('ensembl_gene_id', 'external_gene_name')#'ensembl_transcript_id', 'go_id', 'external_gene_name', 'description')

  # Run the query
  annot <- getBM(attributes=attributeNames,
                 filters = filterType,
                 values = filterValues,
                 mart = ensembl)
  
  return(annot)
}

removeDuplicates <- function(counts.data) {
  # genes <- gsub("\\..+", "", rownames(counts.data))
  # isDup <- duplicated(genes)
  # dup <- results$genes[isDup]
  # results[results$genes%in%dup,]
  counts.data$gene <- gsub("\\..+", "", rownames(counts.data))
  counts.data <- unique(counts.data)
  toRemove <- grep('PAR_Y', rownames(counts.data))
  if(!is_empty(toRemove)) {
    print(c("to remove: ", toRemove))
    counts.data <- counts.data[-toRemove,]
  }
  rownames(counts.data) <- counts.data$gene
  counts.data <- counts.data[,-ncol(counts.data)]
  return(counts.data)
}

# Remove lowly expressed genes
remove.lowly.expressed.genes <- function(counts.data) {
  ## remove lowly expressed genes
  ## here genes that have 4 counts or less in at least 5% of patients
  nb.patients <- as.integer(dim(counts.data)[2])
  thresh <- counts.data > 4
  # rowSums(thresh)
  keep <- rowSums(thresh) >= floor(0.05 * nb.patients)
  # summary(keep)
  counts.keep <- counts.data[keep, ]
  # logcounts.keep <- log(counts.keep)
  
  return(counts.keep)
}

# Create dge object (an EdgeR object) from a RNA-seq matrix
create.dge.object <- function(counts, features, annot) {
  dgeObj <- DGEList(counts[-1,], samples = features, genes = annot[,2])
  return(dgeObj)
}

# Compute log-CPM from a dge object
log.cpm.distrib <- function(dgeObj) {
  return(cpm(dgeObj, log = TRUE))
}

# Turns a vector of character into a list that can be displayed in a selectInput object
make.features.list <- function(df) {
  chlst <- c("", colnames(df))
  names(chlst) <- c("", colnames(df))
  return(chlst)
}

# Deduce LRT model from dge object
analyze.de <- function(dgeObj, feature, condition) {
  dgeObj <- estimateCommonDisp(dgeObj)
  dgeObj <- estimateGLMTrendedDisp(dgeObj)
  dgeObj <- estimateTagwiseDisp(dgeObj)
  # plotBCV(dgeObj)
  feature.vector <- dgeObj$samples[, as.character(feature)]
  feature.vector <- feature.vector == condition
  design <- model.matrix(~ feature.vector)
  fit <- glmFit(dgeObj, design)
  lrt <- glmLRT(fit, coef = 2)
  return(lrt)
}

# Find the enriched functions for the de genes
getEnrichedFunctions <- function(df) {
  genes <- df$FDR < 0.01
  names(genes) <- rownames(def)
  pwf <- nullp(genes, "hg19", "ensGene")
  go.results <- goseq(pwf, "hg19","ensGene", test.cats = c("GO:BP"))
  enriched.GO <- go.results$category[p.adjust(go.results$over_represented_pvalue,method="BH")<.01]
  enrichedGoResults <- as.data.frame(merge(go.results, as.data.frame(enriched.GO), by.x="category", by.y="enriched.GO"))
  enrichedFunctions() %>%
    top_n(10, wt=-over_represented_pvalue) %>%
    mutate(hitsPerc=numDEInCat*100/numInCat) %>% test

  return(test)
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


# ============================================================
# UI
# ============================================================
ui <- fluidPage(
  # Title
  # ============================================================
  titlePanel("RNA-seq data Analyzer"),
  
  # Layout (sidebar layout)
  # ============================================================
  sidebarLayout(
    
    # Sidebar
    # --------------------------------------
    sidebarPanel(
    # Subpanel1: Input the data
      tags$h4(style="text-align: center", "Load data"),
      inputPanel(
        ## Input for the counts table
        fileInput(inputId = "counts.file",
                  label = "Select a RNA-seq count table \n "),
        # tags$param("Samples must be columns and genes must be rows. The first line must be a header with the sample names. The first columns must contain the gene identifiers."),
        
        ## Input for the features table
        fileInput(inputId = "features.file",
                  label = "Select the associated clinical features"),
        
        ## checkbox for early annotation
        checkboxInput(inputId = "annotateGenes", 
                      label = "annotated genes (may take some time)",
                      value = FALSE),
        
        ## Action button
        # verbatimTextOutput(outputId="file"),
        actionButton(inputId = "explore", label = "Exploratory Analysis")
    ),
    
    # Subanel2: Choose the feature and condition to test against the other group(s)
    br(),
    tags$h4(style="text-align: center", "Differential Expression"),
    wellPanel(
      ## Clinical feature for patient groups
      selectInput(
        inputId = "mdsGroupingFeature",
        label = "Choose a clinical feature to group your samples",
        choices = ""
      ),
      
      ## Radio button to choose condition to test
      conditionalPanel(
        condition = "input.mdsGroupingFeature != ''",
        radioButtons(
          inputId = "DE.condition",
          label = "Choose a condition to test against the others",
          choices = ""
        )
      ),
      
      ## Action buttons
        actionButton(inputId = "run.de", label = "Differential Expression Analysis"),
        actionButton(inputId = "save.de", label = "Save model"),
        ## Input for the counts table
        fileInput(inputId = "loadDE",
                label = "Load model (rds file)")
    )
  ), 
  
  # Main Panel
  # --------------------------------------
  mainPanel(
    # Subpanel1: Explore the data
    wellPanel(
      "Data Exploration",
      ## First row: simple size and logCPM distributions
      fluidRow(column(6,
                      plotOutput(outputId = "lib.sizes")), # , width="90%")),
               column(
                 6,
                 plotOutput(outputId = "norm.logcpm.distrib") # , width="90%")
               )),
      br(),
      ## Second row: MDS and most variable genes heatmap
      fluidRow(
        column(6,
               align = "center",
               fluidRow(plotOutput(outputId = "mds")),
               fluidRow(
                 conditionalPanel(
                   condition = "input.explore != 0",
                   actionButton(inputId = "run.glimma", label = "Glimma MDS")
                 )
               )),
        column(6,
               plotOutput(outputId = "most.variable.genes"))
      ),
      htmlOutput("glimma") # , inline=TRUE)
    ),
    
    # Subpanel2: Differential expression
    wellPanel(
      "Differential Gene Expression",
      fluidRow(verbatimTextOutput(outputId = "de.summary")),
      br(),
      fluidRow(column(6,
                      plotOutput(outputId = "de.pval")),
               column(6,
                      plotOutput(outputId = "de.fdr"))),
      br(),
      fluidRow(column(
        6,
        fluidRow(plotOutput(outputId = "de.smear")),
        fluidRow(actionButton(inputId = "runGlimmaSmear", label = "Glimma Smear"))
      ),
      column(6,
             fluidRow(plotOutput(outputId = "de.volcano")),
             fluidRow(actionButton(inputId = "runGlimmaVolcano", label = "Glimma Volcano Plot"))
      ),
      br(),
      ),
      br(),
      # fluidRow(
      #   plotOutput(outputId = "deEnrichedFunctions")
      # ),
      # htmlOutput(outputId="de.volcano")),
      htmlOutput("glimmaSmear"),
      htmlOutput("glimmaVolcano"),
      tableOutput(outputId = "enrichedFunctions")#,
      # fluidRow(plotOutput(outputId = "deEnrichedFunctions"))

    ))
))
                

# ==============================
# SERVER
# ==============================
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 1000 * 1024 ^ 2)
  ## Get some nice colours using Rcolorbrewer
  rdybu.palette <- brewer.pal(11, "RdYlBu")
  set2.palette <- brewer.pal(8, "Set2")
  mds.colors <- colorRampPalette(set2.palette)
  heatmap.colors <-
    colorRampPalette(rdybu.palette)
  
  # Reactive Values
  # ============================================================
  # The properties of the input count table
  file <- reactive(input$features.file)
  
  # The counts table
  counts <- eventReactive(input$explore, {
    data.matrix(remove.lowly.expressed.genes(removeDuplicates(getRawData(input$counts.file$datapath))))
  })
  
  # The features table
  features <- eventReactive(input$explore, {
    as.data.frame(getRawData(input$features.file$datapath))
  })
  
  # The list of clinical features
  chlst <- eventReactive(input$explore, {
    make.features.list(as.data.frame(features()))
  })
  
  annotatedGenes <- eventReactive(input$explore, {
    if(input$annotateGenes) {
      print("Sending Query to Biomart")
      annotateGenes(counts())
    }
    else {
      NULL
    }
  })
  
  # keep <- eventReactive(input$explore, {
  #   intersect(names(counts()), names(t(features())))
  # })
  
  # The dge obect created by EdgeR after loading the data
  # and removing lowly expressed genes
  dgeObj <- eventReactive(input$explore, {
    create.dge.object(as.matrix(counts()), features(), annotatedGenes())
  })
  
  dgeObj.norm <- eventReactive(dgeObj(), {
    calcNormFactors(dgeObj())
  })
  
  # The logCPM distribution of reads deduced from the dge object
  ## choose a clinical feature to group your data
  # output$features <- renderPrint(chlst())
  observeEvent(chlst(), {
    updateSelectInput(session, "mdsGroupingFeature", choices = chlst())
  })
  
  mdsGroupingFeature <-
    reactive(input$mdsGroupingFeature)
  
  logCPM <- eventReactive(dgeObj(), {
    cpm(dgeObj(), log = TRUE)
  })
  
  norm.logCPM <- eventReactive(dgeObj.norm(), {
    cpm(dgeObj.norm(), log = TRUE)
  })
  # limma.groups <- reactive({
  #   groups <- as.factor(merge(dge()$sample, features(),by="row.names")[,as.character(input$mdsGroupingFeature)])
  # })
  #
  observeEvent(input$run.glimma, {
    if (input$mdsGroupingFeature != "") {
      glMDSPlot(
        dgeObj.norm(),
        labels = rownames(dgeObj.norm()$samples),
        groups = dgeObj.norm()$samples[, as.character(input$mdsGroupingFeature)]
      )
    } else {
      glMDSPlot(dgeObj.norm(), labels = rownames(dgeObj.norm()$samples))
    }
  })
  
  lrt <- eventReactive(input$run.de, {
    analyze.de(dgeObj.norm(),
               mdsGroupingFeature(),
               as.character(input$DE.condition))
  })
  
  de.df <- eventReactive(lrt(), {
    as.data.frame(topTags(lrt(), n = Inf))
  })
  
  deDecideTest <- eventReactive(de.df(), {
    decideTestsDGE(lrt(), p = 0.05, adjust = "BH")
  })
  
  deChangedGenes <-
    eventReactive(deDecideTest(), {
      rownames(dgeObj.norm()[as.logical(deDecideTest()), ])
    })

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


  observeEvent(input$save.de, {
    lrt <- lrt()
    saveRDS(lrt, file="./model.rds")
  })
  
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
  
  ## size of libraries
  output$lib.sizes <- renderPlot({
    barplot(
      dgeObj()$samples$lib.size,
      names = colnames(dgeObj()),
      las = 2,
      main = "Barplot of library sizes"
    )
  })
  
  ## log-CPM distribution of libraries
  output$norm.logcpm.distrib <- renderPlot({
    boxplot(
      norm.logCPM(),
      xlab = "",
      ylab = "Log2 counts per million",
      las = 2,
      main = "Boxplots of logCPMs (normalised)"
    )
    abline(h = median(norm.logCPM()), col = "blue")
  })
  
  ## MDS
  ## compute and plot the MDS
  # output$selected.mds.group <- renderPrint({
  #   as.numeric(as.factor(features()[,as.character(input$mdsGroupingFeature)]))
  # })
  output$mds <- renderPlot({
    if (input$mdsGroupingFeature != "") {
      # plotMDS(dge(), labels=labels, groups=group, folder="mds")
      feat <-
        as.factor(dgeObj.norm()$samples[, as.character(input$mdsGroupingFeature)])
      plotMDS(dgeObj.norm(),
              pch = 16,
              main = "MDS",
              col = set2.palette[as.numeric(feat)])
      legend(
        "topleft",
        col = set2.palette,
        pch = 16,
        legend = levels(as.factor(feat))
      )
    } else {
      plotMDS(dgeObj.norm(), pch = 16, main = "MDS")
    }
  })
  
  output$glimma <- renderUI({
    includeHTML("glimma-plots/MDS-Plot.html")
  })
  
  # Plot the most variable genes
  output$most.variable.genes <- renderPlot({
    ## get the 500 most variable genes and their name
    var_genes <- apply(logCPM(), 1, var)
    select_var <-
      names(sort(var_genes, decreasing = TRUE))[1:500]
    highly_variable_lcpm <-
      logCPM()[select_var, ]
    ## plot them
    heatmap.2(
      highly_variable_lcpm,
      col = rev(heatmap.colors(50)),
      trace = "none",
      main = "Top 500 most variable genes across samples",
      # ColSideColors=col.cell,
      scale = "row"
    )
  })
  
  #-- Panel Differential Expression ---------------------------------------
  
  ## update the possible conditions to choose for DE
  # observeEvent(mdsGroupingFeature(), {
  #   updateRadioButtons(session, "DE.control", choices=unique(dgeObj()$samples[,as.character(input$mdsGroupingFeature)]))
  # })
  observeEvent(mdsGroupingFeature(), {
    updateRadioButtons(session, "DE.condition", choices = unique(dgeObj()$samples[, as.character(input$mdsGroupingFeature)]))
  })
  
  ## display the summary of the DE analysis
  output$de.summary <- renderPrint({
    summary(deDecideTest())
  })
  
  # Plot the histogram of p-values for DE genes
  output$de.pval <- renderPlot({
    # print(head(de.df()))
    hist(de.df()$PValue,
         las = 2,
         main = "Histogram of p-values",
         xlab = "p-values")
  })
  
  ## plot the histogram of FDR
  output$de.fdr <- renderPlot({
    hist(de.df()$FDR,
         las = 2,
         main = "Histogram of FDR",
         xlab = "FDR")
    abline(v = 0.05, col = "blue")
    
  })
  
  # Plot the smear with selected genes
  output$de.smear <- renderPlot({
    plotSmear(lrt(), de.tags = deChangedGenes(), main = "Smear plot")
    # points(de.changed.genes()$logCPM, de.changed.genes()$logFC, col="red")#, labels = significant$SYMBOL,col="red")
  })
  observeEvent(input$runGlimmaSmear, {
    output$glimmaSmear <- renderUI({
      glMDPlot(lrt(),
               counts=counts()[-1,],
               xlab="logFC", 
               ylab="-log(FDR)", 
               status=deDecideTest(), 
               anno=annotatedGenes(), 
               groups=as.factor(dgeObj.norm()$samples[, as.character(input$mdsGroupingFeature)]),
               side.main="external_gene_nama"
               )
    })
  })

  # Plot the volcano plot with selected genes
  output$de.volcano <- renderPlot({
    signif <- -log10(de.df()$FDR)
    plot(
      de.df()$logFC,
      signif,
      pch = 16,
      main = "Volcano plot",
      xlab = "logFC",
      ylab = expression("-log"["10"] * "(FDR)")
    )
    points(de.df()[deChangedGenes(), "logFC"], -log10(de.df()[deChangedGenes(), "FDR"]), pch = 16, col = "red")
    abline(v = -2, col = "blue")
    abline(v = 2, col = "blue")
    abline(h = 1.3, col = "green")
  })
  
  observeEvent(input$runGlimmaVolcano, {
    print(head(de.df()))
    print(head(dgeObj.norm()$counts))
    output$glimmaVolcano <- renderUI({
      ga2=data.frame(GeneID=rownames(de.df()), rownames=rownames(de.df()))
      glXYPlot(x= de.df()$logFC,
               y=-log10(de.df()$FDR),
               xlab="logFC", 
               ylab="-log(FDR)", 
               status=as.numeric(de.df()$FDR <= 0.05), 
               anno=ga2)
    })
    
    # glXYPlot(df$logFC,-log10(df$FDR),
    #          xlab="logFC",
    #          ylab="-log(FDR)",
    #          status=as.numeric(df$FDR <= 0.05),
    #          anno=ga2)
  })
  output$annotatedGenes <- renderTable(head(de.df()))

#   output$deEnrichedFunctions <- renderPlot({
#     list(ggplot(enrichedFunctions(), aes_string(x=hitsPerc,
#                  y=term,
#                  colour=over_represented_pvalue,
#                  size=numDEInCat) 
#            +geom_point()
#            + expand_limits(x=0)
#            + labs(x="Hits (%)", y="GO term", colour="p value", size="Count")
#     ))
#   })
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