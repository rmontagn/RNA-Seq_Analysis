# ==============================
# Load packages
# ==============================
list.of.packages <-
  c("ggplot2", "Rcpp", "rmarkdown", "gplots", "RColorBrewer")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
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
new.packages <-
  list.of.bioconductor.packages[!(list.of.bioconductor.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  BiocManager::install(new.packages)

rm(list = ls())

library(edgeR)
# library(DESeq)
# library(limma)
library(Glimma)
library(gplots)
# library(org.Mm.eg.db)
library(RColorBrewer)
# library(Homo.sapiens)
library(shiny)
library(GO.db)

# ============================================================
# SERVER FUNCTIONS
# ============================================================

# Get raw counts table
get.raw.data <- function(counts.file) {
  seqdata <-
    read.delim(
      as.character(counts.file),
      stringsAsFactors = FALSE,
      header = F,
      sep = "\t"
    )
  colnames(seqdata) <- seqdata[1, ]
  rownames(seqdata) <- seqdata[, 1]
  counts.data <- seqdata[-1, -1]
  # counts.data <- data.matrix(counts.data)
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
create.dge.object <- function(counts, features) {
  dgeObj <- DGEList(counts, samples = features)
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

# dt <- as.data.frame(dt)
# results <- merge(results, dt, by=0)
# rownames(results) <- results$Row.names
# # anno <- as.data.frame(cbind(results$Row.names, results$Row.names))
# head(results)
# test <- c("group1", "group2")[group+1]
# glMDPlot(lrt(), counts = dgeObj(), groups = test,
#          xlab="logFC", ylab="-log(FDR)", status=results$survival1, anno=results)


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
        # width = "100%"),
        
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
      
      ## Action button
      actionButton(inputId = "run.de", label = "Differential Expression Analysis")
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
      # htmlOutput(outputId="de.volcano")),
      htmlOutput("glimmaSmear"),
      htmlOutput("glimmaVolcano")
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
    data.matrix(remove.lowly.expressed.genes(get.raw.data(input$counts.file$datapath)))
  })
  
  # The features table
  features <- eventReactive(input$explore, {
    as.data.frame(get.raw.data(input$features.file$datapath))
  })
  
  # The list of clinical features
  chlst <- eventReactive(input$explore, {
    make.features.list(as.data.frame(features()))
  })
  
  # keep <- eventReactive(input$explore, {
  #   intersect(names(counts()), names(t(features())))
  # })
  
  # The dge obect created by EdgeR after loading the data
  # and removing lowly expressed genes
  dgeObj <- eventReactive(input$explore, {
    create.dge.object(as.matrix(counts()), features())
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
  
  de.decide.test <- eventReactive(de.df(), {
    decideTestsDGE(lrt(), p = 0.05, adjust = "BH")
  })
  
  de.changed.genes <-
    eventReactive(de.decide.test(), {
      rownames(dgeObj.norm()[as.logical(de.decide.test()), ])
    })
  
  
  # Pots
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
    summary(de.decide.test())
  })
  
  # Plot the histogram of p-values for DE genes
  output$de.pval <- renderPlot({
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
    plotSmear(lrt(), de.tags = de.changed.genes(), main = "Smear plot")
    # points(de.changed.genes()$logCPM, de.changed.genes()$logFC, col="red")#, labels = significant$SYMBOL,col="red")
  })
  observeEvent(input$runGlimmaSmear, {
    output$glimmaSmear <- renderUI({
      glMDPlot(lrt(), xlab="logFC", ylab="-log(FDR)", status=de.decide.test())
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
    points(de.df()[de.changed.genes(), "logFC"], -log10(de.df()[de.changed.genes(), "FDR"]), pch = 16, col = "red")
    abline(v = -2, col = "blue")
    abline(v = 2, col = "blue")
    abline(h = 1.3, col = "green")
  })
  observeEvent(input$runGlimmaVolcano, {
    output$glimmaVolcano <- renderUI({
      df <- de.df()[rownames(de.decide.test()@.Data),]
      signif <- -log10(df$FDR)
      glXYPlot(df$logFC, y=signif, xlab="logFC", ylab="-log(FDR)", status=de.decide.test())
    })

  })

}

shinyApp(ui = ui, server = server)
