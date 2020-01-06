# RNA-Seq_Analysis
A shiny application for RNA-seq data anlysis

## Provide inputs

![interface01](readme/interface01.png)

Load a count table and a feature table:  
- The count table must be a tab-separated file. The rows must be genes and the columns samples. The first row must correspond to the table header with the sample names. The first column must correspond to the genes names or keys from an ontology.  
- The feature table must be a tab-separated file. The rows must be the samples and the columns a feature characterizing the dataset. The first row must correspond to the name of each feature. The first column must correspond to the sample names. The sample names must be exactly the same as in the count table.  
The sample in the count table and the feature table must be exactly the same.

## Data exploration

Click `Exploratory Analysis`. The interface displays the sample sizes, logCPM, the MDS plot and a heatmap of the 500 most variable genes. 

In the section `Choose a clinical feature to group your samples`, select a clinical feature. This will classify the data and the MDS will be interactively colored according to the different groups. If you wish more functionalities, click `Glimma plot`. This will open in a web browser an interactive MDS plot. The package used to generate this MDS, `glimma`, is different from the one used in the interface, `EdgeR`, so the MDS can be a little different.

## Differential gene expression

Once a feature has been chosen to classify the samples (cf previous section), the list of its possible values are displayed in the section `Choose a condition to test against the others`. You can choose one of these values. Click `Differential Expression Analysis`. The interface will compute the differential gene expression between the group of samples corresponding to this value and the other ones. Other possibilities need to be implemented.