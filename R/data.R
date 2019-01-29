##############################################
# The deveroper and the maintainer:          #
#   Lang Ho Lee (lhlee@bwh.harvard.edu)      #
#   Sasha A. Singh (sasingh@bwh.harvard.edu) #
##############################################

#' Randomly generated example datasets for XINA users.
#' A dataset containing the XINA clustering results.
#'
#' \itemize{
#'   \item aligned. XINA clustering results aligned by conditions
#'   \item data_column. Column names for data matrix
#'   \item out_dir. Not available in this example dataset
#'   \item nClusters. The number of user-desired clusters.
#'   It's 30 in the example.
#'   \item max_cluster. The number of clusters found in the dataset.
#'   It's 21 in the example.
#'   \item chosen_model. The chosen covariance model
#'   for the example dataset.  It's VEI in the example
#'   \item optimal_BIC. BIC at the optimized clustering.
#'   It's 29473.57 in the example
#'   \item condition. The experimental conditions in the dataset.
#'   \item color_for_condition. The default color
#'   for the conditions that will be used in XINA plot drawing.
#'   \item color_for_clusters. The default color
#'   for the clusters that will be used in XINA clustering plot.
#'   \item norm_method. The used normalization method
#'   to standardize the input data. It's "sum_normalization"
#'   in the example.
#' }
#'
#' @format A list with the example XINA clustering result
#' @name example_clusters
NULL

#' Protein-protein interaction resource downloaded from HPRD DB
#' A data frame containing HRPD protein-protein interaction data
#'
#' \itemize{
#'   \item gene_symbol_1. Gene name interacting with gene name
#'   in 'gene_symbol_2'
#'   \item gene_symbol_2. Gene name interacting with gene name
#'   in 'gene_symbol_1'
#'   \item Experiment_type. Experimental or computational methods
#'   supporting the interaction
#' }
#'
#' @format A data frame containing HRPD protein-protein interaction data
#' @source \url{http://www.hprd.org/}
#' @name hprd_ppi
NULL

#' Protein-protein interaction resource downloaded from STRING DB
#' for XINA's example dataset
#' A data frame containing protein-protein interactions
#'
#' \itemize{
#'   \item gene_symbol_1. Gene name interacting with gene name
#'   in 'gene_symbol_2'
#'   \item gene_symbol_2. Gene name interacting with gene name
#'   in 'gene_symbol_1'
#'   \item PPI_Source. Data original source
#' }
#'
#' @format A data frame containing STRING protein-protein interaction data
#' @source \url{https://string-db.org/}
#' @name string_example
NULL

#' Previously processed xina analysis using XINA's random example data
#' A list containing 'xina_analysis' results
#'
#' \itemize{
#'   \item All_network. PPI network of all the input proteins
#'   \item Sub_network. A list containing PPI networks of each clusters
#'   \item Data. XINA clustering results. See \link[XINA]{xina_clustering}
#'   \item Nodes. A list of proteins in each cluster
#'   \item Conditions. A list of experimental condition of proteins
#'   in each cluster
#'   \item Titles.  A list of plot titles for XINA plotting
#'   \item out_dir. A directory path storing XINA network analysis results
#'   \item is_stringdb. False = different PPI DB and TRUE = STRING DB
#' }
#'
#' @format A data frame containing STRING protein-protein interaction data
#' @source XINA
#' @name xina_result_example
NULL


#' A character vector containing 19,396 human gene descriptions
#' This is for the randome data generation of XINA
#'
#' \itemize{
#'   \item Human gene description corresponding to 'gn' vector
#' }
#'
#' @format A character vector containing 19,396 human gene descriptions
#' @source \url{https://www.ncbi.nlm.nih.gov/gene}
#' @name gn_desc
NULL

#' A character vector containing 19,396 human genes
#' This is for the randome data generation of XINA
#'
#' \itemize{
#'   \item Characters of human genes
#' }
#'
#' @format A character vector containing 19,396 human genes
#' @source \url{https://www.ncbi.nlm.nih.gov/gene}
#' @name gn
NULL
