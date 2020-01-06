# ==============================
# Load packages
# ==============================
list.of.packages <- c("ggplot2", "Rcpp", "rmarkdown", "gplots", "RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

list.of.bioconductor.packages <- c("limma", "Glimma", "edgeR", "DESeq", "org.Hs.eg.db", "Homo.sapiens", "GO.db")
new.packages <- list.of.bioconductor.packages[!(list.of.bioconductor.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

rm(list=ls())

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

# ==============================
# SERVER FUNCTIONS
# ==============================
# Get the raw counts table
get.raw.data <- function(counts.file) {
  seqdata <- read.delim(as.character(counts.file), stringsAsFactors = FALSE, header=F, sep="\t")
  colnames(seqdata) <- seqdata[1,]
  rownames(seqdata) <- seqdata[,1]
  counts.data <- seqdata[-1, -1]
  # counts.data <- data.matrix(counts.data)
  return(counts.data)
}

# Remove lowly expressed genes
remove.lowly.expressed.genes <- function(counts.data) {
  # Remove lowly expressed genes
  # Here genes that have 4 counts or less in at least 5% of patients 
  nb.patients = as.integer(dim(counts.data)[2])
  thresh <- counts.data > 4
  # rowSums(thresh)
  keep <- rowSums(thresh) >= floor(0.05*nb.patients)
  # summary(keep)
  counts.keep <- counts.data[keep,]
  # logcounts.keep <- log(counts.keep)
  
  return(counts.keep)
}

# Create dge object (an EdgeR object) from a RNA-seq matrix
create.dge.object <- function(counts, features) {
  dgeObj <- DGEList(counts, samples=features)
  return(dgeObj)
}

# Compute the log-CPM from a dge object
log.cpm.distrib <- function(dgeObj) {
  return(cpm(dgeObj,log=TRUE))
}

# Turns a vector of character into a list that can be displayed in a selectInput object
make.features.list <- function(df) {
  chlst <- c("", colnames(df))
  names(chlst) <- c("", colnames(df))
  return(chlst)
}

analyze.de <- function(dgeObj, feature, condition) {
  dgeObj <- estimateCommonDisp(dgeObj)
  dgeObj <- estimateGLMTrendedDisp(dgeObj)
  dgeObj <- estimateTagwiseDisp(dgeObj)
  # plotBCV(dgeObj)
  feature.vector <- dgeObj$samples[,as.character(feature)]
  feature.vector <- feature.vector == condition
  design <- model.matrix(~ feature.vector)
  fit <- glmFit(dgeObj, design)
  lrt <- glmLRT(fit, coef=2)
  return(lrt)
}

# ==============================
# UI
# ==============================
ui <- fluidPage(
  titlePanel("RNA-seq data Analyzer"),
  sidebarLayout(
    sidebarPanel(
      
      # Panel1: Input the data
      inputPanel(
        #Input for the counts table
        fileInput(
         inputId = "counts.file",
         label = "Select a RNA-seq count table",
          # width = "25%"
        ),
        #tags$param("Samples must be columns and genes must be rows. The first line must be a header with the sample names. The first columns must contain the gene identifiers."),
        
        #Input for the features table  
        fileInput(
          inputId = "features.file",
          label = "Select the associated clinical features",
          # width = "100%"
        ),
        
        # Action button
        # verbatimTextOutput(outputId="file"),
        actionButton(inputId="explore", label="Exploratory Analysis")
      ),

      # Panel2: choose the grouping feature for MDS and subsequent DE analysis
      wellPanel(
        selectInput(inputId="mds.grouping.feature", 
                    label="Choose a clinical feature to group your samples",
                    choices = ""
        ),
        actionButton(inputId="run.glimma", label="Glimma plot"),
      ),
      
      # Panel3: Choose the condition to test against the other group(s)
      wellPanel(
        radioButtons(inputId="DE.condition", label="Choose a condition to test against the others", choices=""),
        actionButton(inputId="run.de", label="Differential Expression Analysis")
      )
    
    ),
    
    mainPanel(
      wellPanel(
        "Data Exploration",
        fluidRow(
          column(6,
            plotOutput(outputId="lib.sizes")#, width="90%")
            ),
          column(6,
                 plotOutput(outputId="norm.logcpm.distrib")#, width="90%")
          )
        ),
        br(),
        fluidRow(
          column(6,
                  plotOutput(outputId="mds")
          ),
          column(6,
                  plotOutput(outputId="most.variable.genes")
          )
        ),
        htmlOutput("glimma", inline=TRUE)
      ),

      wellPanel(
        "Differential Gene Expression",
        fluidRow(
          verbatimTextOutput(outputId="de.summary")
        ),
        br(),
        fluidRow(
          column(6,
                 plotOutput(outputId="de.pval"),
          ),
          column(6,
                 plotOutput(outputId="de.fdr")
          )
        ),
        br(),
        fluidRow(
          column(6,
                 plotOutput(outputId="de.smear")
          ),
          column(6,
                 plotOutput(outputId="de.volcano")
          )
        )
      )
    )
  )
)

# ==============================
# SERVER
# ==============================
server <- function(input, output, session){
  options(shiny.maxRequestSize=1000*1024^2)
  ## Get some nice colours using Rcolorbrewer
  rdybu.palette <- brewer.pal(11,"RdYlBu")
  set2.palette <- brewer.pal(8, "Set2")
  mds.colors <- colorRampPalette(set2.palette)
  heatmap.colors <- colorRampPalette(rdybu.palette)
  
  # ===============
  # Get the reactive values
  # ===============
  # The properties of the input count table
  file <- reactive(input$features.file)
  
  # The counts table
  counts <- eventReactive(input$explore, {
    data.matrix(
      remove.lowly.expressed.genes(
        get.raw.data(input$counts.file$datapath)
      )
    )
  })
  
  # The features table
  features <- eventReactive(input$explore, {
    as.data.frame(
          get.raw.data(input$features.file$datapath)
    )
  })
  
  # The list of clinical features
  chlst <- eventReactive(input$explore, {
    make.features.list(
      as.data.frame(
        features()
      )
    )
  })
  
  # keep <- eventReactive(input$explore, {
  #   intersect(names(counts()), names(t(features())))
  # })
  
  # The dge obect created by EdgeR after loading the data
  # and removing lowly expressed genes
  dgeObj <- eventReactive(input$explore, { 
    create.dge.object(
      as.matrix(counts()), features()
    )
  })
  
  dgeObj.norm <- eventReactive(dgeObj(), { 
    calcNormFactors(dgeObj())
  })
  
  # The logCPM distribution of reads deduced from the dge object  
  mds.grouping.feature <- reactive(input$mds.grouping.feature)
  
  logCPM <- eventReactive(dgeObj(), {
    cpm(dgeObj(),log=TRUE)
  })
  
  norm.logCPM <- eventReactive(dgeObj.norm(), {
    cpm(dgeObj.norm(), log=TRUE)
  })
  # limma.groups <- reactive({
  #   groups <- as.factor(merge(dge()$sample, features(),by="row.names")[,as.character(input$mds.grouping.feature)])
  # })
  # 
  glimma <- eventReactive(input$run.glimma, {
    glMDSPlot(dgeObj.norm(), labels=rownames(dgeObj.norm()$samples), groups=dgeObj.norm()$samples[,as.character(input$mds.grouping.feature)])
  })
  
  lrt <- eventReactive(input$run.de, {
    analyze.de(dgeObj.norm(), mds.grouping.feature(), as.character(input$DE.condition))
  })
  
  de.df <- eventReactive(lrt(), {
    as.data.frame(topTags(lrt(), n=Inf))
  })
  
  de.decide.test <- eventReactive(de.df(), {
    decideTestsDGE(lrt(), p=0.05, adjust="BH")
  })

  de.changed.genes <- eventReactive(de.decide.test(), {
    rownames(dgeObj.norm()[as.logical(de.decide.test()),])
  })
  
  # ===============
  # Plot the results
  # ===============
  # Plot the counts table file propertiess
  output$file <- renderPrint(file())  

  # Plot the log-CPM distribution of libraries
  output$lib.sizes <- renderPlot({
    barplot(dgeObj()$samples$lib.size, names=colnames(dgeObj()), las=2, main="Barplot of library sizes")
  })
  
  output$norm.logcpm.distrib <- renderPlot({
    boxplot(norm.logCPM(), xlab="", ylab="Log2 counts per million", las=2, main="Boxplots of logCPMs (normalised)")
    abline(h=median(norm.logCPM()), col="blue")
  })

  # Plot the MDS of the dataset
  ## choose a clinical feature to group your data
  # output$features <- renderPrint(chlst())
  observeEvent(chlst(), {
    updateSelectInput(session, "mds.grouping.feature", choices=chlst())
  })

  ## compute and plot the MDS
  # output$selected.mds.group <- renderPrint({
  #   as.numeric(as.factor(features()[,as.character(input$mds.grouping.feature)]))
  # })
  output$mds <- renderPlot({
    if(input$mds.grouping.feature != ""){
      # plotMDS(dge(), labels=labels, groups=group, folder="mds")
      feat <- as.factor(dgeObj.norm()$samples[,as.character(input$mds.grouping.feature)])
      plotMDS(dgeObj.norm(), pch=16, main="MDS", col=set2.palette[as.numeric(feat)])
      legend("topleft", col=set2.palette, pch=16, legend=levels(as.factor(feat)))} else {
        plotMDS(dgeObj.norm(), pch=16, main="MDS")
        }
  })
  output$glimma <- renderUI({
    glimma()
  })

  # Plot the most variable genes
  output$most.variable.genes <- renderPlot({
    ## get the 500 most variable genes and their name
    var_genes <- apply(logCPM(), 1, var)
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
    highly_variable_lcpm <- logCPM()[select_var,]
    ## plot them
    heatmap.2(highly_variable_lcpm, 
              col=rev(heatmap.colors(50)),
              trace="none", 
              main="Top 500 most variable genes across samples",
              # ColSideColors=col.cell,
              scale="row"
    )
  })

  # Plot the results of differentiall genes expression
  
  ## update the possible conditions to choose for DE
  # observeEvent(mds.grouping.feature(), {
  #   updateRadioButtons(session, "DE.control", choices=unique(dgeObj()$samples[,as.character(input$mds.grouping.feature)]))
  # })
  observeEvent(mds.grouping.feature(), {
    updateRadioButtons(session, "DE.condition", choices=unique(dgeObj()$samples[,as.character(input$mds.grouping.feature)]))
  })
  
  ## display the summary of the DE analysis
  output$de.summary <- renderPrint({
    summary(de.decide.test())
  })

  # Plot the histogram of p-values for DE genes
  output$de.pval <- renderPlot({
    hist(de.df()$PValue, las=2, main="Histogram of p-values", xlab="p-values")
  })
  
  ## plot the histogram of FDR
  output$de.fdr <- renderPlot({
    hist(de.df()$FDR, las=2, main="Histogram of FDR", xlab="FDR")
    abline(v=0.05, col="blue")
  })
  
  # Plot the smear with selected genes
  output$de.smear <- renderPlot({
    plotSmear(lrt(), de.tags=de.changed.genes(), main="Smear plot")
    # points(de.changed.genes()$logCPM, de.changed.genes()$logFC, col="red")#, labels = significant$SYMBOL,col="red")
  })
  
  # Plot the volcano plot with selected genes
  output$de.volcano <- renderPlot({
    signif <- -log10(de.df()$FDR)
    plot(de.df()$logFC, signif, pch=16, main="Volcano plot", xlab="logFC", ylab=expression("-log"["10"]*"(FDR)"))
    points(de.df()[de.changed.genes(),"logFC"],-log10(de.df()[de.changed.genes(),"FDR"]),pch=16,col="red")
    abline(v=-2, col="blue")
    abline(v=2, col="blue")
    abline(h=1.3, col="green")
  })

}

shinyApp(ui=ui, server=server)