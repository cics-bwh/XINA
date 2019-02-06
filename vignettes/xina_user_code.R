## ----setup, include=FALSE---------------------------------------------------------------------------------------------
options(width=120)
knitr::opts_chunk$set(
    echo=TRUE
)

## ----installation, eval=FALSE, warning=FALSE--------------------------------------------------------------------------
#  # Install from Bioconductor
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("XINA")
#  
#  # Install from Github
#  install.packages("devtools")
#  library(devtools)
#  install_github("langholee/XINA")

## ----import libraries, warning=FALSE----------------------------------------------------------------------------------
library(XINA)

## ----import required packages, eval=FALSE, warning=FALSE--------------------------------------------------------------
#  install.packages("igraph")
#  install.packages("ggplot2")
#  BiocManager::install("STRINGdb")

## ----example random dataset, warning=FALSE----------------------------------------------------------------------------
# Generate random multiplexed time-series data
random_data_info <- make_random_xina_data()

# The number of proteins
random_data_info$size

# Time points
random_data_info$time_points

# Three conditions
random_data_info$conditions

## ----check randomly generated data files, warning=FALSE---------------------------------------------------------------
Control <- read.csv("Control.csv", check.names=FALSE, stringsAsFactors = FALSE)
Stimulus1 <- read.csv("Stimulus1.csv", check.names=FALSE, stringsAsFactors = FALSE)
Stimulus2 <- read.csv("Stimulus2.csv", check.names=FALSE, stringsAsFactors = FALSE)

head(Control)
head(Stimulus1)
head(Stimulus2)

## ----data matrix, warning=FALSE---------------------------------------------------------------------------------------
head(Control[random_data_info$time_points])

## ----fix random seed, warning=FALSE-----------------------------------------------------------------------------------
set.seed(0)

## ----set up for the clustering, warning=FALSE-------------------------------------------------------------------------
# Data files
data_files <- paste(random_data_info$conditions, ".csv", sep='')
data_files

# time points of the data matrix
data_column <- random_data_info$time_points
data_column

## ----XINA clustering, warning=FALSE-----------------------------------------------------------------------------------
# Run the model-based clusteirng
clustering_result <- xina_clustering(data_files, data_column=data_column, nClusters=20)

## ----XINA with k-means clustering, warning=FALSE----------------------------------------------------------------------
clustering_result_km <- xina_clustering(data_files, data_column=data_column, nClusters=20, chosen_model='kmeans')

## ----XINA clustering plot, warning=FALSE------------------------------------------------------------------------------
library(ggplot2)
theme1 <- theme(title=element_text(size=8, face='bold'),
                axis.text.x = element_text(size=7),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank())
plot_clusters(clustering_result, ggplot_theme=theme1)

## ----XINA condition composition, warning=FALSE------------------------------------------------------------------------
theme2 <- theme(legend.key.size = unit(0.3, "cm"),
                legend.text=element_text(size=5),
                title=element_text(size=7, face='bold'))
condition_composition <- plot_condition_compositions(clustering_result, ggplot_theme=theme2)
tail(condition_composition)

## ----XINA comigration search, warning=FALSE---------------------------------------------------------------------------
classes <- as.vector(clustering_result$condition)
classes

all_cor <- alluvial_enriched(clustering_result, classes)
head(all_cor)

## ----XINA comigration search with comigration size filter, eval=TRUE, warning=FALSE-----------------------------------
cor_bigger_than_5 <- alluvial_enriched(clustering_result, classes, comigration_size=5)
head(cor_bigger_than_5)

## ----STRING DB set up, eval=FALSE, warning=FALSE----------------------------------------------------------------------
#  library(STRINGdb)
#  string_db <- STRINGdb$new( version="10", species=9606, score_threshold=0, input_directory="" )
#  string_db
#  
#  xina_result <- xina_analysis(clustering_result, string_db)

## ----plotting PPI networks of all the clusters, eval=FALSE, warning=FALSE---------------------------------------------
#  # XINA network plots labeled gene names
#  xina_plot_all(xina_result, clustering_result)

