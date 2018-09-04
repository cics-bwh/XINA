#############################################
# The deveroper and the maintainer:         #
#   Lang Ho Lee (lhlee@bwh.harvard.edu)     #
#   Sasha A. Singh (sasingh@bwh.harvard.edu)#
#############################################

# Mute warnings
options(warn=1)

#' @title find_similar_clusters
#' @description Compare clusters and find similar ones
#' @param clustering_result A list containing XINA clustering results.
#' See \link[XINA]{xina_clustering}
#' @param threshold Pearson's r threshold to find similar ones
#' @return Write a csv file containing similar clustering information
#' based on the given Pearson's R threshold
#' @import utils
#'
find_similar_clusters <- function(clustering_result, threshold=0.95)
{
  Cluster <- NULL
  clustering_dir <- clustering_result$out_dir
  nClusters <- clustering_result$nClusters
  data_column <- clustering_result$data_column
  column_numbers <- seq_len(length(data_column))
  super_ds <- clustering_result$clusters

  # Compare clusters using PCC
  df_mean_vectors <- data.frame()
  for (i in seq_len(nClusters))
  {
    cluster_data <- subset(super_ds, Cluster==i)
    if (nrow(cluster_data) != 0) {
      mean_vector <- c()
      for (j in seq_len(length(data_column)))
      {
        mean_vector <- c(mean_vector, mean(cluster_data[data_column[j]][[1]]))
      }
      tmp <- cbind(mean_vector)
      colnames(tmp) <- i
      rownames(tmp) <- data_column
      if (nrow(df_mean_vectors)==0){
        df_mean_vectors <- rbind(df_mean_vectors, tmp)
      } else {
        df_mean_vectors <- cbind(df_mean_vectors, tmp)
      }
    }
  }
  cor_matrix <- stats::cor(df_mean_vectors)
  df_similar_cluster <- data.frame()
  for (i in seq_len(nrow(cor_matrix)))
  {
    for (j in i:ncol(cor_matrix))
    {
      if (i != j)
      {
        if (cor_matrix[i,j] >= threshold){
          tmp <- rbind(c(i,j,cor_matrix[i,j]))
          colnames(tmp) <- c("C1","C2","PCC")
          df_similar_cluster <- rbind(df_similar_cluster, tmp)
        }
      }
    }
  }
  write.csv(df_similar_cluster,paste(clustering_dir,"/Similar_Clusters.csv", sep=''),row.names = FALSE)
}

#' @title extract_data_column
#' @description Extract data column names from XINA clustering result
#' @param col_head_of_clustering Column names of XINA clustering result
#' @return A vector containing column names of  data matrix
extract_data_column <- function(col_head_of_clustering){
  cluster_info_head <-c("Accession", "Description", "Condition",
                        "Key", "Cluster")
  data_column <- c()
  for (col_head in col_head_of_clustering){
    if (is.na(match(col_head, cluster_info_head))){
      data_column <- c(data_column, col_head)
    }
  }
  return(data_column)
}

#' @title load_previous_results
#' @description Get previous XINA clustering results to R space
#' @param clustering_dir The directory path of XINA clustering results
#' @param data_column A vector containing column names of  data matrix
#' @param fp_clusters File path of XINA clustering results
#' @return Comma-separated file containing aligned XINA clustering results.
#' @import utils
#' @export
#' @examples
#' # Load XINA's example data
#' data(xina_example)
#' write.csv(example_clusters$aligned,"xina_clusters_aligned.csv")
#' write.csv(example_clusters$clusters,"xina_clusters.csv")
#'
#' # Reload the clustering result
#' example_clusters_reloaded <- load_previous_results(".")
#'
load_previous_results <- function(clustering_dir=getwd(), data_column=NULL,
                                  fp_clusters="xina_clusters.csv") {
  # Retrieve clustering results
  file_clusters <- paste(clustering_dir, "/", fp_clusters, sep='')
  df_clusters <- read.csv(file_clusters, check.names=FALSE,
                          stringsAsFactors = FALSE)
  df_aligned <- organize_clusters(clustering_dir, df_clusters, file_out=FALSE)
  if (is.null(data_column)) {
    data_column <- extract_data_column(colnames(df_clusters))
  }
  conditions <- unique(df_clusters$Condition)
  condition_colors <- get_color_for_nodes()[seq_len(length(conditions))]
  names(condition_colors) <- conditions
  nClusters <- max(df_clusters$Cluster)
  color_for_clusters <- get_colors(nClusters, set='original')
  names(color_for_clusters) <- seq_len(nClusters)
  # Recall previous results
  mIMT_result <- list()
  mIMT_result$clusters <- df_clusters
  mIMT_result$aligned <- df_aligned
  rownames(mIMT_result$clusters) <- mIMT_result$clusters$Key
  mIMT_result$data_column <- data_column
  mIMT_result$out_dir <- clustering_dir
  mIMT_result$nClusters <- nClusters
  mIMT_result$max_cluster <- nClusters
  mIMT_result$chosen_model <- NA
  mIMT_result$optimal_BIC <- NA
  mIMT_result$condition <- conditions
  mIMT_result$color_for_condition <- condition_colors
  mIMT_result$color_for_clusters <- color_for_clusters
  return(mIMT_result)
}


#' @title organize_clusters
#' @description Organize XINA clustering information by gene name
#' @param clustering_dir The directory path of XINA clustering results
#' @param super_ds XINA clusters
#' @param file_out If it is TRUE, it writes the aligned clustering informaion
#' to "xina_clusters_aligned.csv" file.
#' @import utils
#' @return Comma-separated file containing aligned XINA clustering results.
organize_clusters <- function(clustering_dir=getwd(), super_ds, file_out=TRUE)
{
  Accession <- Condition <- NULL
  conditions <- unique(super_ds$Condition)
  genes <- unique(super_ds$Accession)

  df_cluster <- data.frame()
  for (i in seq_len(length(genes)))
  {
    ds_gene <- subset(super_ds, Accession == genes[i])
    clusters <- c()
    for (j in seq_len(length(conditions)))
    {
      cluster <- subset(ds_gene, Condition==conditions[j])
      if (nrow(cluster)>0)
      {
        cluster_number <- as.integer(cluster$Cluster)
      } else {
        cluster_number <- 0
      }
      clusters <- c(clusters, cluster_number)
    }
    tmp <- c(as.character(ds_gene$Accession[1]),
             as.character(ds_gene$Description[1]), clusters)
    df_cluster <- rbind(df_cluster, data.frame(t(tmp),
                                               stringsAsFactors = FALSE))
  }
  colnames(df_cluster) <- c("Gene name", "Description", as.vector(conditions))
  if (file_out) {
    write.csv(df_cluster,"xina_clusters_aligned.csv",row.names = FALSE)
  }
  return(df_cluster)
}

#' @title xina_clustering
#' @description Clustering multiplexed time-series omics data to
#' find co-abundance profiles
#' @param f_names A vector containing input file (.csv) paths
#' @param data_column A vector containing column names
#' (1st row of the input file) of data matrix
#' @param out_dir A directory path for saving clustering results.
#' (default: out_dir=getwd())
#' @param nClusters The number of desired maximum clusters
#' @param chosen_model You can choose a specific model rather than
#' testing all the models that are available in mclust.
#' \link[mclust]{mclustModelNames}
#' If you want k-means clustering instead of the model-based clustering,
#' use "kmeans" here.
#' @param norm Default is "sum_normalization".
#' Sum-normalization is to divide the data matrix by row sum.
#' If you want to know more about sum-normalization,
#' see https://www.ncbi.nlm.nih.gov/pubmed/19861354.
#' "zscore" is to calculate Z score for each protein.
#' See \link[base]{scale}.
#' @return a plot containing a BIC plot in current working directory
#' and a list containing below information:
#'      \tabular{rl}{
#'       \strong{Item} \tab \strong{Description}\cr
#'       clusters \tab XINA clustering results\cr
#'       aligned \tab XINA clustering results aligned by ID\cr
#'       data_column \tab Data matrix column names\cr
#'       out_dir \tab The directory path containing XINA results\cr
#'       nClusters \tab The number of clusters desired by user\cr
#'       max_cluster \tab The number of clusters optimized by BIC\cr
#'       chosen_model \tab The used covariance model for model-based clustering\cr
#'       optimal_BIC \tab BIC of the optimized covariance model\cr
#'       condition \tab Experimental conditions of the user input data\cr
#'       color_for_condition \tab Colors assigned to each experimental conditions
#'       which is used for condition composition plot\cr
#'       color_for_clusters \tab Colors assigned to each clusters
#'       which is used for XINA clustering plot\cr
#'       norm_method \tab Used normalization method\cr
#'      }
#' @import mclust
#' @import utils
#' @export
#' @examples
#' # Generate random multiplexed time-series data
#' random_data_info <- make_random_xina_data()
#'
#' # Data files
#' data_files <- paste(random_data_info$conditions, ".csv", sep='')
#'
#' # time points of the data matrix
#' data_column <- random_data_info$time_points
#'
#' # mclust requires the fixed random seed to get reproduce the clustering results
#' set.seed(0)
#'
#' # Run the model-based clustering to find co-abundance profiles
#' example_clusters <- xina_clustering(data_files, data_column=data_column,
#' nClusters=30)
#'
#' # Run k-means clustering to find co-abundance profiles
#' example_clusters <- xina_clustering(data_files, data_column=data_column,
#' nClusters=30,
#' chosen_model="kmeans")
#'
xina_clustering <- function(f_names, data_column, out_dir=getwd(), nClusters=20,
                            norm="sum_normalization", chosen_model="")
{
  Cluster <- NULL
  clustering_dir <- getwd()
  # Merge data
  superset <- generate_superset(f_names, data_column, norm=norm)
  rownames(superset) <- superset$Key
  # Run Mclust
  if (chosen_model == "kmeans") {
    #superset <- na.omit(superset[data_column]) # listwise deletion of missing
    # mydata <- scale(mydata) # standardize variables
    fit <- stats::kmeans(superset[data_column], nClusters)
    superset$Cluster <- fit$cluster
    highest_bic <- NA
  } else {
    if (chosen_model == "") {
      mc_Sds <- Mclust(superset[data_column],G=2:nClusters)
      chosen_model <- mc_Sds$modelName
    } else {
      mc_Sds <- Mclust(superset[data_column],G=2:nClusters, modelNames=chosen_model)
    }
    superset$Cluster <- mc_Sds$class
    plot.Mclust(mc_Sds, what="BIC")
    # Write messages from Mclust
    sink(paste("./mclust_summary.txt", sep=''))
    summSds<-summary(mc_Sds,parameters=TRUE)
    print(summSds)
    sink()
    highest_bic <- mc_Sds$bic
  }
  max_cluster <- max(superset$Cluster)

  # Remove empty clusters
  superset_nr <- data.frame()
  cluster_num <- 0
  for (i in seq_len(max_cluster)){
    sub_cluster <- subset(superset, Cluster==i)
    if (nrow(sub_cluster) > 0){
      cluster_num <- cluster_num + 1
      sub_cluster[,"Cluster"] <- cluster_num
      superset_nr <- rbind(superset_nr, sub_cluster)
    }
  }
  write.csv(superset_nr, file="xina_clusters.csv", row.names = FALSE)
  # Collect cluster information by proteins
  aligned_clusters <- organize_clusters(clustering_dir, superset_nr)
  conditions <- unique(superset_nr$Condition)
  max_cl <- max(superset_nr$Cluster)
  # Set default condition colors
  condition_colors <- get_color_for_nodes()[seq_len(length(conditions))]
  names(condition_colors) <- conditions
  # Set default cluster colors
  color_for_clusters <- get_colors(nClusters, set='original')
  names(color_for_clusters) <- seq_len(nClusters)
  # return the return list
  return(list("clusters"=superset_nr,
              "aligned"=aligned_clusters,
              "data_column"=data_column,
              "out_dir"=clustering_dir,
              "nClusters"=nClusters,
              "max_cluster"=max_cluster,
              "chosen_model"=chosen_model,
              "optimal_BIC"=highest_bic,
              "condition"=conditions,
              "color_for_condition"=condition_colors,
              "color_for_clusters"=color_for_clusters,
              "norm_method"=norm))
}

#' @title generate_superset
#' @description Merge input kinetics files
#' @param f_names A vector of .csv file paths containing kinetics data
#' @param data_column A vector of column names containing data matrix
#' @param delim The delimiter of input file (default is ',')
#' @param norm The normalization method. It should be one of
#' c('sum_normalization', 'zscore').
#' Default is 'sum_normalization'.
#' @import tools
#' @import utils
#' @return A data frame containing kinetics data obtained from files
#' in the f_names vector
generate_superset <- function(f_names, data_column, delim=",",
                              norm="sum_normalization")
{
  SuperDS <- data.frame()
  for (i in seq_len(length(f_names)))
  {
    f_name <- f_names[i]
    condition <- file_path_sans_ext(f_name)
    df_csv <- read.csv(file=f_name, check.names=FALSE, sep=delim)
    if (is.na(match("Accession",colnames(df_csv)))){
      stop(paste("Your input file named '",f_name,
                 "' doesn't have 'Accession' column"))
    }
    if (is.na(match("Description",colnames(df_csv)))){
      stop(paste("Your input file named '",f_name,
                 "' doesn't have 'Description' column"))
    }
    # remove zero-sum rows
    df_filtered <- subset(df_csv, rowSums(df_csv[data_column])>0)
    # remove na rows
    df_filtered <- df_filtered[stats::complete.cases(df_filtered[data_column]), ]
    # normalize the data matrix
    if (norm=="sum_normalization") {
      # Sum-normalization
      df_normalized <- cbind(df_filtered[data_column] /
                               as.matrix(rowSums(df_filtered[data_column])),
                             Accession=df_filtered$Accession,
                             Description=df_filtered$Description, Condition=condition)
    } else if (norm=="zscore") {
      # Z-score
      df_normalized <- data.frame()
      for (k in seq_len(nrow(df_filtered))) {
        tmp <- as.vector(scale(as.numeric(df_filtered[k, data_column])))
        names(tmp) <- data_column
        df_normalized <- rbind(df_normalized, tmp)
        colnames(df_normalized) <- data_column
      }
      df_normalized["Accession"] <- df_filtered$Accession
      df_normalized["Description"] <- df_filtered$Description
      df_normalized["Condition"] <- condition
    } else {
      stop("normalization method should be one of these c(sum_normalization, zscore)")
    }
    # Add unique keys
    df_normalized["Key"] <- paste(df_normalized$Accession, "@", condition, sep="")
    # Remove duplicated keys
    df_normalized=df_normalized[!duplicated(df_normalized$Key),]
    # append the normalized row
    SuperDS <-rbind(SuperDS, df_normalized)
  }
  return(SuperDS)
}
