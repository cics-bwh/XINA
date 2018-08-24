#############################################
# The deveroper and the maintainer:         #
#   Lang Ho Lee (lhlee@bwh.harvard.edu)     #
#   Sasha A. Singh (sasingh@bwh.harvard.edu)#
#############################################

# Mute warnings
options(warn=1)

#' @title make_random_xina_data
#' @description Generate random proteomics dataset for testing XINA
#' 'make_random_xina_data' will make random proteomics data for XINA test.
#' The generated data will have three conditions and seven time points,
#' c("0hr", "2hr", "6hr", "12hr", "24hr", "48hr", "72hr").
#' @param n The number of proteins for one condition. Default is 500.
#' @param mtor If it is TRUE (default), mTOR pathway genes will be significant.
#' If it is FALSE, randomly selected genes will be significant
#' in first three conditions.
#' @param time_points A vector containing time points of the data matrix
#' @param conditions A vector containing condition information,
#' for example normal, disease and drug treated disase.
#' @return Three comma-separated files containing time-series data for XINA
#' @import utils
#' @export
#' @examples
#' make_random_xina_data()
#' g1 <- read.csv("Control.csv", check.names=FALSE,
#' stringsAsFactors = FALSE)
#' g2 <- read.csv("Stimulus1.csv", check.names=FALSE,
#' stringsAsFactors = FALSE)
#' g3 <- read.csv("Stimulus2.csv", check.names=FALSE,
#' stringsAsFactors = FALSE)
#'
#' head(g1)
#' head(g2)
#' head(g3)
#'
make_random_xina_data <- function(n=500, mtor=TRUE,
                                  time_points=c("0hr", "2hr", "6hr", "12hr",
                                                "24hr", "48hr", "72hr"),
                                  conditions=c("Control", "Stimulus1",
                                               "Stimulus2")) {
  # gene_info <- gn_desc <- gn <- NULL
  # Get gene names and descriptions
  data("gene_info", envir=environment())
  gn <- gn
  gn_desc <- gn_desc
  # Get mTOR pathway data
  if (mtor) {
    list_intended_genes <- get_mTOR_proteins(time_points, conditions)
    nrow_intended <- nrow(list_intended_genes[[1]])
  } else {
    nrow_intended <- 0
  }
  # Calculate required numbers
  num_total <- n - nrow_intended
  num_shared <- round(num_total*0.8)
  num_unique <- num_total - num_shared
  limit_gn <- length(gn) - nrow_intended
  # Get random data
  if (mtor) {
    list_data <- get_random_data(time_points, conditions, num_shared+num_unique,
                                 percent.sign=0)
  } else {
    list_data <- get_random_data(time_points, conditions, num_shared+num_unique)
  }
  # Get random gene names
  shared_idx <- sample(x=seq_len(limit_gn), size=num_shared)
  for (i in seq_len(length(conditions))) {
    selected_idx <- shared_idx
    for (j in seq_len(num_unique)) {
      tmp_idx <- 0
      while (TRUE) {
        tmp_idx <- sample(x=seq_len(limit_gn), size=1)
        if (is.na(match(tmp_idx, selected_idx))) {
          break
        }
      }
      selected_idx <- c(selected_idx, tmp_idx)
    }
    if (mtor) {
      data_out <- rbind(cbind(Accession=gn[selected_idx],
                              Description=gn_desc[selected_idx],
                              list_data[[i]]),
                        list_intended_genes[[i]])
    } else {
      data_out <- cbind(Accession=gn[selected_idx],
                        Description=gn_desc[selected_idx], list_data[[i]])
    }
    write.csv(data_out, paste(conditions[i],".csv",sep=''), row.names = FALSE)
  }
  return(list(size=n, time_points=time_points, conditions=conditions))
}

#' @title get_random_data
#' @description Get randomized time-series data
#' @param time_points A vector containing time points of the data matrix
#' @param conditions A vector containing condition information,
#' for example normal, disease and drug treated disase.
#' @param num_total The number of total proteins to be generated
#' @param percent.sign Percentage of differentially expressed proteins.
#' Ignored when equal=FALSE.
#' @param equal If equal is TRUE, all the conditions will have numbers
#' between 0 and 1.
#' If it is  FALSE, the first three conditions will have different ranges.
#' First condition will have numbers from 0.3 to 0.4.
#' Second condition will have numbers from 0.6 to 0.8.
#' Third condition will have numbers from 0.3 to 0.5.
#' Other conditions will have numbers from 0 to 1.
#' @return A list containing ramdomly generated data matrix
#'
get_random_data <- function(time_points, conditions, num_total,
                            percent.sign=0.1, equal=TRUE){
  if (equal) {
    significant_idx <- sample(x=seq_len(num_total),
                              size=round(num_total*percent.sign))
  } else {
    significant_idx <- seq_len(num_total)
  }
  list_data <- list()
  for (i in seq_len(num_total)) {
    for (j in seq_len(length(conditions))) {
      if (is.na(match(2,significant_idx))) {
        expression_profile <- c(0.5, stats::runif(length(time_points)-1))
      } else {
        if (j==1) {
          expression_profile <-
            c(0.5, stats::runif(length(time_points)-1,0.3,0.4))
        } else if (j==2) {
          expression_profile <-
            c(0.5, stats::runif(length(time_points)-1,0.6,0.8))
        } else if (j==3) {
          expression_profile <-
            c(0.5, stats::runif(length(time_points)-1,0.3,0.5))
        } else {
          expression_profile <-
            c(0.5, stats::runif(length(time_points)-1))
        }
      }
      names(expression_profile) <- time_points
      if (length(list_data)<j) {
        list_data[[j]] <- round(expression_profile,4)
      } else {
        tmp <- rbind(list_data[[j]], round(expression_profile,4))
        colnames(tmp) <- time_points
        rownames(tmp) <- NULL
        list_data[[j]] <- tmp
      }
    }
  }
  return(list_data)
}

#' @title get_mTOR_proteins
#' @description Get mTOR pathway genes
#' @param time_points A vector containing time points of the data matrix
#' @param conditions A vector containing condition information,
#' for example normal, disease and drug treated disase.
#' @import utils
#' @return A vector containing mTOR pathway gene names
get_mTOR_proteins <- function(time_points, conditions) {
  # gene_info <- gn_desc <- gn <- NULL
  data("gene_info", envir=environment())
  gn <- gn
  gn_desc <- gn_desc
  mtor_gn <- c("AKT1","AKT1S1","BRAF","CAB39","CAB39L","CBX2",
               "DDIT4","EIF4B","EIF4E","EIF4E2","EIF4EBP1",
               "HIF1A","IGF1","IKBKB","IRS1","MAPK1","MAPK3","MLST8",
               "MTOR","PDPK1","PIK3CA","PIK3CB","PIK3CD","PIK3CG",
               "PIK3R1","PIK3R2","PIK3R3","PIK3R5","PRKAA1","PRKAA2",
               "PRKCA","PRKCB","PRKCG","PTEN","RHEB","RICTOR","RPS6",
               "RPS6KA1","RPS6KA2","RPS6KA3","RPS6KA6","RPS6KB1",
               "RPS6KB2", "RPTOR","RRAGA","RRAGB","RRAGC","RRAGD","RYBP",
               "STK11","STRADA", "TNF","TSC1","TSC2","ULK1","ULK2","VEGFA")
  mtor_desc <- gn_desc[match(mtor_gn, gn)]
  list_data <- get_random_data(time_points, conditions, length(mtor_gn),
                               equal=FALSE)
  list_return <- list()
  for (i in seq_len(length(conditions))) {
    list_return[[i]] <- cbind(Accession=mtor_gn, Description=mtor_desc,
                              list_data[[i]])
  }
  return(list_return)
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


#' @title get_unknown_ppi_nodes
#' @description Get proteins with no known interactions within the cluster based on the used protein-protein interaction database source
#' @param xina_result A list containing XINA network analysis results. See \link[XINA]{xina_analysis}
#' @param cl the clustering number of XINA clustering results. See \link[XINA]{xina_clustering}
#' @return A data frame containing proteins with no known interactions within the cluster based on the used protein-protein interaction database source
#' @import igraph
#' @export
#' @examples
#' # load XINA example data
#' data(xina_example)
#'
#' # load the previously processed XINA analysis results
#' # if you want to learn how to run 'xina_analysis', please see \link[XINA]{xina_analysis}
#' data(xina_result_example)
#'
#' # Extract unknown PPI nodes in the cluster #1
#' get_unknown_ppi_nodes(xina_result_example, 1)
#'
get_unknown_ppi_nodes <- function(xina_result, cl){
  subnet <- xina_result$Sub_network[[cl]]
  isolated_nodes <- V(subnet)[degree(subnet)==0]
  df_subnet <- data.frame(Name=isolated_nodes$name, Condition=isolated_nodes$condition, Color=isolated_nodes$vertex.color)
  return(df_subnet)
}

#' @title add_legend
#' @description Add plot legend and locate it outside of a network plot
#' @param legend_location Network centrality score matrix
#' @param ... Numeric, complex, or logical vectors.
#' @return a legend to a plot
#'
add_legend <- function(legend_location="bottomright", ...) {
  opar <- graphics::par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
                        mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(graphics::par(opar))
  graphics::plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  graphics::legend(legend_location, ...)
}

#' @title xina_plot_single
#' @description xina_plot_single draws protein-protein interaction network plot for given 'protein_list'.
#'
#' @param xina_result A list containing XINA network analysis results. See \link[XINA]{xina_analysis}
#' @param protein_list A vector of gene names to draw a protein-protein interaction network graph.
#' @param centrality_type 'centrality_type' should be one of
#' c('Degree', 'Eigenvector', 'Hub', 'Authority', 'Closeness', 'Betweenness')
#'      \tabular{rl}{
#'       \strong{Centrality score} \tab \strong{igraph function}\cr
#'       Degree \tab \link[igraph]{degree}\cr
#'       Eigenvector \tab \link[igraph]{eigen_centrality}\cr
#'       Hub \tab \link[igraph]{hub_score}\cr
#'       Authority \tab \link[igraph]{authority_score}\cr
#'       Closeness \tab \link[igraph]{closeness}\cr
#'       Betweenness \tab \link[igraph]{betweenness}\cr
#'      }
#' @param layout_specified This can change network layout.
#' 'layout_specified' should be one of c('sphere', 'star', 'gem', 'tree', 'circle', 'random', 'nicely').
#' XINA's layouts are based on igraph's layout. See \link[igraph]{layout_}
#'      \tabular{rl}{
#'       \strong{Layout} \tab \strong{igraph layout name}\cr
#'       sphere \tab \link[igraph]{layout_on_sphere}\cr
#'       star \tab \link[igraph]{layout_as_star}\cr
#'       gem \tab \link[igraph]{layout_with_gem}\cr
#'       tree \tab \link[igraph]{layout_as_tree}\cr
#'       circle \tab \link[igraph]{layout_in_circle}\cr
#'       random \tab \link[igraph]{layout_randomly}\cr
#'       nicely \tab \link[igraph]{layout_nicely}\cr
#'      }
#' Default is 'layout_nicely' of igraph
#' @param vertex_label_flag
#' If vertex_label_flag is TRUE (default), igraph network graphs will be labeled by gene names
#' If vertex_label_flag is FALSE, igraph network graphs will be drawn without labels
#' @param main Title of network figure.  IF it is NULL (default), it will be the number of plotted proteins
#' @param vertex.label.color Color of labels. Default is black
#' @param vertex.color Color of nodes. Default is pink.
#' @param edge.color Color of edges. Default is pink.
#' @param vertex.label.dist Distance between node and label. Default is 0.6
#' @param vertex.label.cex Size of labels  Default is 0.8
#' @param edge.arrow.size Size of edges  Default is 0.4
#' @param vertex.size Size of nodes  Default is 10
#' @param vertex.shape You can choose node shape. Default is 'sphere'.  See \link[igraph]{shapes}
#' @param legend_location If centrality_type is chosen,
#' 'xina_plot_single' adds the color legend guiding rank of nodes based on the centrality score.
#' Default is 'bottomright', but you can choose one of these 'bottomright', 'bottom', 'bottomleft',
#' 'left', 'topleft', 'top', 'topright', 'right' and 'center'.
#' @param num_breaks 'num_breaks' is the number of ranks based on network centrality. Default is 5.
#' @param digits_round_up See \link[base]{Round}
#' @param flag_simplify If it is TRUE (default), XINA will exclude unconnected proteins
#' @param flag_legend If it is TRUE, a legend will be printed out together.
#' @return A PNG file (XINA_Cluster_Networks.png) displaying protein-protein interaction network plots
#' of all the clusters and a list containing XINA network analysis results
#' @import igraph
#' @export
#' @examples
#' ## the following code is to show how it works quickly
#' ## load XINA example data
#' data(xina_example)
#'
#' ## load the previously processed XINA analysis results
#' # if you want to learn how to run 'xina_analysis', please see \link[XINA]{xina_analysis}
#' data(xina_result_example)
#'
#' # get gene names that are clustered to #21 in "Stimulus2" condition
#' subgroup <- subset(example_clusters$aligned, Stimulus2==21)
#' protein_list <- subgroup$`Gene name`
#'
#' # Calculate protein-protein interaction network
#' xina_plot_single(xina_result_example, protein_list)
#'
#' # Calculate protein-protein interaction network and Eigenvector centrality
#' eigen_info <- xina_plot_single(xina_result_example, protein_list, centrality_type='Eigenvector')
#'
xina_plot_single <- function(xina_result, protein_list, centrality_type=NULL,
                             layout_specified='', vertex_label_flag=TRUE, main=NULL,
                             vertex.label.color='black', vertex.color=NA, edge.color='darkgray',
                             vertex.label.dist=.6, vertex.label.cex=0.8,
                             edge.arrow.size=.4, vertex.size=10, vertex.shape='sphere',
                             legend_location='bottom', num_breaks=5, digits_round_up=5,
                             flag_simplify=TRUE, flag_legend=TRUE) {
  if (is.null(protein_list)) {
    stop("XINA needs at least one valid protein to get PPI network")
  }
  if (length(protein_list[!is.na(protein_list)])==0) {
    stop("XINA needs at least one valid protein to get PPI network")
  }
  net_all <- xina_result$All_network
  vertices <- toupper(as.vector(V(net_all)$name))

  protein_list <- toupper(protein_list)
  nodes <- c()
  for (k in seq_len(length(protein_list))) {
    nodes <- c(nodes, match(protein_list[k], vertices))
  }
  # Get a subnetwork using the matched proteins and simplify the network if flag_simplify is TRUE
  subnet <- induced_subgraph(net_all, nodes[!is.na(nodes)])
  if (flag_simplify) {
    isolated_nodes <- V(subnet)[degree(subnet)==0]
    subnet_simplified <- delete.vertices(subnet, isolated_nodes)
    subnet <- subnet_simplified
  }
  if (length(V(subnet))==0) {
    plot_NA()
    warning("No available nodes in the subnetwork. This can be casued by flag_simplify.  Try xina_plot_single again with flag_simplify=FALSE")
  } else {
    # Plot a network
    if (is.null(main)){
      main <- paste(length(nodes[!is.na(nodes)]),"Proteins",centrality_type)
    }
    # Set up the vertex label
    if (vertex_label_flag) {
      vertex_label <- V(subnet)$name
    } else {
      vertex_label <- NA
    }
    # Select network graph layout
    if (layout_specified=='sphere') {
      layout_selected <- layout_on_sphere(subnet)
    } else if (layout_specified=='star') {
      layout_selected <- layout_as_star(subnet)
    } else if (layout_specified=='gem') {
      layout_selected <- layout_with_gem(subnet)
    } else if (layout_specified=='tree') {
      layout_selected <- layout_as_tree(subnet)
    } else if (layout_specified=='circle') {
      layout_selected <- layout_in_circle(subnet)
    } else if (layout_specified=='random') {
      layout_selected <- layout_randomly(subnet)
    } else if (layout_specified=='nicely') {
      layout_selected <- layout_nicely(subnet)
    } else {
      layout_selected <- get_layout(subnet)
    }
    if (is.null(centrality_type)) {
      vector_colors <- vertex.color
    } else {
      centrality_score <- calculate_centrality_scores(subnet, centrality_type)
      centrality_results <- data.frame(Gene_Name=V(subnet)$name,
                                       Description=V(subnet)$desc,
                                       Score=centrality_score,
                                       Rank=rank_centrality(centrality_score,centrality_type,num_breaks=num_breaks))
      tgc <- get_stats(centrality_results)
      rbPal <- colorRampPalette(c("white", "red"))
      vector_colors <- rbPal(num_breaks)[centrality_results$Rank]
    }
    plot(subnet, vertex.label.color="black", layout=layout_selected,
         vertex.label.dist=vertex.label.dist, vertex.label.cex=vertex.label.cex,
         vertex.label=vertex_label, edge.arrow.size=edge.arrow.size,
         vertex.size=vertex.size, vertex.color=vector_colors,
         edge.color=edge.color, vertex.shape=vertex.shape, main=main)
    # return centrality scores
    if (!is.null(centrality_type)) {
      if (flag_legend) {
        # Add a legend and
        add_legend(legend_location=legend_location, legend=paste("<",round(tgc$max,digits_round_up)), pch=20,
                   col=rbPal(length(tgc$max)), horiz=TRUE, bty='n', cex=vertex.label.cex, title=centrality_type)
      }
      return(centrality_results)
    }
  }
}

#' @title xina_plot_bycluster
#' @description xina_plot_bycluster is to draw protein-protein interaction network plots of each cluster
#' @param xina_result A list containing XINA network analysis results. See \link[XINA]{xina_analysis}
#' @param clustering_result A list containing XINA clustering results. See \link[XINA]{xina_clustering}
#' @param cl Cluster number in the XINA clustering results
#' @param condition Default is 'all', which means use all the proteins to draw graphs.
#' If you specify the experimental condition name used for XINA clustering,
#' @param flag_legend If it is TRUE, a legend will be printed out together.
#' @param centrality_type 'centrality_type' should be one of
#' c('Degree', 'Eigenvector', 'Hub', 'Authority', 'Closeness', 'Betweenness')
#'      \tabular{rl}{
#'       \strong{Centrality score} \tab \strong{igraph function}\cr
#'       Degree \tab \link[igraph]{degree}\cr
#'       Eigenvector \tab \link[igraph]{eigen_centrality}\cr
#'       Hub \tab \link[igraph]{hub_score}\cr
#'       Authority \tab \link[igraph]{authority_score}\cr
#'       Closeness \tab \link[igraph]{closeness}\cr
#'       Betweenness \tab \link[igraph]{betweenness}\cr
#'      }
#' @param flag_simplify If it is TRUE (default), XINA will exclude unconnected proteins
#' @param layout_specified This can change network layout.
#' 'layout_specified' should be one of c('sphere', 'star', 'gem', 'tree', 'circle', 'random', 'nicely').
#' XINA's layouts are based on igraph's layout. See \link[igraph]{layout_}
#'      \tabular{rl}{
#'       \strong{Layout} \tab \strong{igraph layout name}\cr
#'       sphere \tab \link[igraph]{layout_on_sphere}\cr
#'       star \tab \link[igraph]{layout_as_star}\cr
#'       gem \tab \link[igraph]{layout_with_gem}\cr
#'       tree \tab \link[igraph]{layout_as_tree}\cr
#'       circle \tab \link[igraph]{layout_in_circle}\cr
#'       random \tab \link[igraph]{layout_randomly}\cr
#'       nicely \tab \link[igraph]{layout_nicely}\cr
#'      }
#' Default is 'layout_nicely' of igraph
#' @param vertex_label_flag
#' If vertex_label_flag is TRUE (default), igraph network graphs will be labeled by gene names
#' If vertex_label_flag is FALSE, igraph network graphs will be drawn without labels
#' @param vertex.color Color of nodes. Default is pink.
#' @param edge.color Color of edges. Default is pink.
#' @param vertex.label.dist Distance between node and label. Default is 0.6
#' @param vertex.label.cex Size of labels  Default is 0.8
#' @param edge.arrow.size Size of edges  Default is 0.4
#' @param vertex.size Size of nodes  Default is 10
#' @param vertex.shape You can choose node shape. Default is 'sphere'.  See \link[igraph]{shapes}
#' @param legend_location If centrality_type is chosen,
#' xina_plot_single add the color legend guiding rank of nodes based on the centrality score.
#' Default is 'bottomright', but you can choose one of these 'bottomright', 'bottom', 'bottomleft',
#' 'left', 'topleft', 'top', 'topright', 'right' and 'center'.
#' @param flag_unknown_only If this is TRUE, 'xina_plot_bycluster' will plot proteins that do not have any protein-protein interaction in the given database
#' @return A PNG file (XINA_Cluster_Networks.png) displaying protein-protein interaction network plots
#' of all the clusters and a list containing XINA network analysis results
#' @return PNG images of PPI network plots of all the clusters
#' @import igraph
#' @export
#' @examples
#' ## the following code is to show how it works quickly
#' ## load XINA example data
#' data(xina_example)
#'
#' ## load the previously processed XINA analysis results
#' # if you want to learn how to run 'xina_analysis', please see \link[XINA]{xina_analysis}
#' data(xina_result_example)
#'
#' # plot cluster #1
#' xina_plot_bycluster(xina_result_example, example_clusters, cl=1)
#'
#' # plot PPI network of Control condition in cluster #1
#' xina_plot_bycluster(xina_result_example, example_clusters, cl=1, condition='Control')
#'
xina_plot_bycluster <- function(xina_result, clustering_result,
                                cl=NULL, condition='all', flag_legend=TRUE,
                                centrality_type=NULL, flag_simplify=TRUE,
                                layout_specified='', vertex_label_flag=TRUE,
                                vertex.label.dist=.6, vertex.label.cex=0.8,
                                edge.arrow.size=.4, vertex.size=10,
                                vertex.shape='sphere', vertex.color='',
                                edge.color='darkgray', legend_location='bottom',
                                flag_unknown_only=FALSE) {
  Condition <- vertex.label.color <- NULL
  # collect clustering information
  nClusters <- clustering_result$nClusters
  data_column <- clustering_result$data_column
  column_numbers <- seq_len(length(data_column))
  uniq_condition <- as.vector(clustering_result$condition)
  max_cluster <- clustering_result$max_cluster
  color_for_nodes <- clustering_result$color_for_condition
  color_for_clusters <- clustering_result$color_for_clusters
  # collect the network information
  titles <- xina_result$Titles
  list_conditions <- xina_result$Conditions
  # get subnetwork
  if (is.null(cl)) {
    stop("You should choose one cluster to draw a network graph")
  }
  subnet <- xina_result$Sub_network[[cl]]
  if (flag_unknown_only) {
    df_subnet <- get_unknown_ppi_nodes(xina_result, cl)
  } else {
    df_subnet <- data.frame(Name=V(subnet)$name, Condition=V(subnet)$condition, Color=V(subnet)$vertex.color)
  }
  if (vertex.color!=''){
    df_subnet$Color <- rep(vertex.color, nrow(df_subnet))
  }
  if (condition =='all') {
    plot_title <- paste("#",cl," (n=", nrow(df_subnet),")", sep='')
  } else {
    if (is.na(match(condition, uniq_condition))) {
      stop("'condition' is not found in XINA clustering results")
    }
    df_subnet <- subset(df_subnet, Condition==condition)
    plot_title <- paste("#",cl," ",condition," (n=", nrow(df_subnet),")", sep='')
  }
  if (nrow(df_subnet)>0) {
    xina_plot_single(xina_result, as.vector(df_subnet$Name), centrality_type=centrality_type,
                     layout_specified=layout_specified, vertex_label_flag=vertex_label_flag,
                     main=plot_title, vertex.label.color=vertex.label.color,
                     vertex.color=df_subnet$Color, edge.color=edge.color,
                     vertex.label.dist=vertex.label.dist, vertex.label.cex=vertex.label.cex,
                     edge.arrow.size=edge.arrow.size, vertex.size=vertex.size,
                     vertex.shape=vertex.shape, legend_location=legend_location,
                     num_breaks=5, flag_simplify=flag_simplify, flag_legend=flag_legend)
  } else {
    plot_NA()
  }
  if (flag_legend) {
    add_legend(legend_location=legend_location, legend=uniq_condition, pch=20,
               col=color_for_nodes, horiz=TRUE, bty='n', cex=vertex.label.cex,
               title="Experimental conditions")
  }
}

#' @title xina_plot_all
#' @description xina_plot_all is to draw protein-protein interaction network plots of all the clusters
#' @param xina_result A list containing XINA network analysis results. See \link[XINA]{xina_analysis}
#' @param clustering_result A list containing XINA clustering results. See \link[XINA]{xina_clustering}
#' @param condition Default is 'all', which means use all the proteins to draw graphs.
#' If you specify the experimental condition name used for XINA clustering,
#' xina_plot_all will draw graphs using specific condition's proteins.
#' @param centrality_type 'centrality_type' should be one of
#' c('Degree', 'Eigenvector', 'Hub', 'Authority', 'Closeness', 'Betweenness')
#'      \tabular{rl}{
#'       \strong{Centrality score} \tab \strong{igraph function}\cr
#'       Degree \tab \link[igraph]{degree}\cr
#'       Eigenvector \tab \link[igraph]{eigen_centrality}\cr
#'       Hub \tab \link[igraph]{hub_score}\cr
#'       Authority \tab \link[igraph]{authority_score}\cr
#'       Closeness \tab \link[igraph]{closeness}\cr
#'       Betweenness \tab \link[igraph]{betweenness}\cr
#'      }
#' @param flag_simplify If it is TRUE (default), XINA will exclude unconnected proteins
#' @param num_breaks 'num_breaks' is the number of ranks based on network centrality. Default is 5.
#' @param layout_specified This can change network layout.
#' 'layout_specified' should be one of c('sphere', 'star', 'gem', 'tree', 'circle', 'random', 'nicely').
#' XINA's layouts are based on igraph's layout. See \link[igraph]{layout_}
#'      \tabular{rl}{
#'       \strong{Layout} \tab \strong{igraph layout name}\cr
#'       sphere \tab \link[igraph]{layout_on_sphere}\cr
#'       star \tab \link[igraph]{layout_as_star}\cr
#'       gem \tab \link[igraph]{layout_with_gem}\cr
#'       tree \tab \link[igraph]{layout_as_tree}\cr
#'       circle \tab \link[igraph]{layout_in_circle}\cr
#'       random \tab \link[igraph]{layout_randomly}\cr
#'       nicely \tab \link[igraph]{layout_nicely}\cr
#'      }
#' Default is 'layout_nicely' of igraph
#' @param vertex_label_flag
#' If vertex_label_flag is TRUE (default), igraph network graphs will be labeled by gene names
#' If vertex_label_flag is FALSE, igraph network graphs will be drawn without labels
#' @param vertex.label.color Color of labels. Default is black
#' @param vertex.color Color of nodes. Default is pink.
#' @param edge.color Color of edges. Default is pink.
#' @param vertex.label.dist Distance between node and label. Default is 0.6
#' @param vertex.label.cex Size of labels  Default is 0.8
#' @param edge.arrow.size Size of edges  Default is 0.4
#' @param vertex.size Size of nodes  Default is 10
#' @param vertex.shape You can choose node shape. Default is 'sphere'.  See \link[igraph]{shapes}
#' @param legend_location If centrality_type is chosen,
#' xina_plot_single add the color legend guiding rank of nodes based on the centrality score.
#' Default is 'bottomright', but you can choose one of these 'bottomright', 'bottom', 'bottomleft',
#' 'left', 'topleft', 'top', 'topright', 'right' and 'center'.
#' @param num_clusters_in_row The number of clusters in a row on the XINA network plot. Default is 5.
#' @param flag_unknown_only If this is TRUE, 'xina_plot_all' will plot proteins that do not have any protein-protein interaction in the given database
#' @param img_size Set the image size. For width=1000 and height=1500,
#' it is img_size=c(1000,1500). Default is c(3000,3000)
#' @param img_qual Set the image resolution. Default is 300.
#' @return PNG images of PPI network plots of all the clusters
#' @import grDevices
#' @import igraph
#' @import graphics
#' @export
#' @examples
#' ## the following code is to show how it works quickly
#' ## load XINA example data
#' data(xina_example)
#'
#' ## load the previously processed XINA analysis results
#' # if you want to learn how to run 'xina_analysis', please see \link[XINA]{xina_analysis}
#' data(xina_result_example)
#'
#' # XINA network plots
#' xina_plot_all(xina_result_example, example_clusters)
#'
#' # XINA network plots for Control condition
#' xina_plot_all(xina_result_example, example_clusters, condition='Control')
#'
xina_plot_all <- function(xina_result, clustering_result, condition='all',
                          centrality_type=NULL, flag_simplify=TRUE, num_breaks=5,
                          layout_specified='', vertex_label_flag=FALSE,
                          vertex.label.color='black', vertex.color='', edge.color=NULL,
                          vertex.label.dist=.6, vertex.label.cex=0.8,
                          edge.arrow.size=.4, vertex.size=10, vertex.shape='sphere',
                          legend_location='bottom', num_clusters_in_row=5,
                          flag_unknown_only=FALSE,
                          img_size=NULL, img_qual=300) {
  # collect clustering information
  nClusters <- clustering_result$nClusters
  uniq_condition <- as.vector(clustering_result$condition)
  max_cluster <- clustering_result$max_cluster
  # image size
  num_clusters_in_row <- xina_result$num_clusters_in_row
  num_clusters_in_col <- as.integer(max_cluster/num_clusters_in_row) + if(max_cluster%%num_clusters_in_row>0){1}else{0}
  if (!is.vector(img_size)){
    unit_size <- 750
    img_size <- c(unit_size*num_clusters_in_row,unit_size*num_clusters_in_col)
  }
  # start png image file
  if (is.null(centrality_type)) {
    category <- condition
  } else {
    category <- paste(condition,"_",centrality_type,sep='')
  }
  out <- paste(xina_result$out_dir,"/","XINA_Cluster_Networks_",category,".png",sep="")
  png(filename=out, width=img_size[1], height=img_size[2], res=img_qual)
  par(mfrow = c(num_clusters_in_col, num_clusters_in_row))
  par(mar=c(1,1,1,1))
  # Draw individual XINA network plots
  for (i in seq_len(max_cluster)){
    if (is.null(edge.color)){
      edge_color <- unique(E(xina_result$Sub_network[[i]])$edge.color)
    } else {
      edge_color <- edge.color
    }
    xina_plot_bycluster(xina_result, clustering_result, cl=i, condition=condition, flag_legend=FALSE,
                        centrality_type=centrality_type, flag_simplify=flag_simplify,
                        edge.color=edge_color, vertex.color=vertex.color,
                        layout_specified=layout_specified, vertex_label_flag=vertex_label_flag,
                        vertex.label.dist=vertex.label.dist, vertex.label.cex=vertex.label.cex,
                        edge.arrow.size=edge.arrow.size, vertex.size=vertex.size,
                        vertex.shape=vertex.shape, legend_location=legend_location,
                        flag_unknown_only=flag_unknown_only)
  }
  dev.off()
}

#' @title plot_enrichment_results
#' @description Plot GO and KEGG enrichment results
#'
#' @param enriched_results GO or KEGG enrichment results.
#' See \link[XINA]{xina_enrichment} and \link[XINA]{xina_enrichment}
#' @param term_description Description of terms to be drawn on Y axis.
#' Default is "term_description" of XINA enrichment results.
#' @param sig_score significant score to plot on X axis. Default is "pvalue".
#' @param num_terms The number of terms to be plotted.  Default is 0, which menas no limit.
#' @param get_log If this is TRUE, 'plot_enrichment_results' will take -log10 of p-values.
#' @return ggplot bar graph
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#' library(STRINGdb)
#'
#' # load XINA example data
#' data(xina_example)
#'
#' # Get STRING database for protein-protein intereaction information
#' string_db <- STRINGdb$new( version="10", species=9606,
#' score_threshold=0, input_directory="" )
#' string_db
#'
#' # XINA analysis with STRING DB
#' xina_result <- xina_analysis(example_clusters, string_db)
#'
#' # Select proteins that showed cluster #1 in the Stimulus2 condition
#' subgroup <- subset(example_clusters$aligned, Stimulus2==1)
#' protein_list <- as.vector(subgroup$`Gene name`)
#'
#' # Enrichment test and get significantly enriched functional terms
#' # that have adjuseted p-value less than 0.1
#' kegg_enriched <- xina_enrichment(string_db, protein_list,
#' enrichment_type = "KEGG", pval_threshold=0.1)
#' plot_enrichment_results(kegg_enriched$KEGG, num_terms=10)
#' }
#'
plot_enrichment_results <- function(enriched_results, term_description="term_description",
                                    sig_score="pvalue", num_terms=0, get_log=TRUE){
  Enriched_score <- NULL
  # requireNamespace("dplyr", quietly = TRUE)
  enriched_results <- enriched_results[!is.na(enriched_results$term_description), ]
  if (get_log) {
    enriched_results["Enriched_score"] <- -log10(enriched_results[sig_score])
    y_lab <- paste("-Log10(",sig_score,")",sep='')
  } else {
    enriched_results["Enriched_score"] <- enriched_results[sig_score]
    y_lab <- paste(sig_score)
  }
  if (term_description!="term_description") {
    enriched_results["term_description"] <- enriched_results[term_description]
  }
  # Limit the results by given the number of terms to display
  if (num_terms > 0 & nrow(enriched_results) >= num_terms) {
    # Order by Enriched_score
    enriched_results <- enriched_results %>% arrange(Enriched_score)
    enriched_results <- enriched_results[(nrow(enriched_results)-num_terms):nrow(enriched_results),]
  }
  enriched_results %>%
    arrange(Enriched_score) %>%
    mutate(term_description = factor(term_description, term_description)) %>%   # reset factor
    ggplot(aes(term_description, Enriched_score)) + geom_bar(stat="identity") + coord_flip() + ylab(y_lab)
}

#' @title default_size
#' @description Calculate image size based on the number of clusters
#' @param max_cluster the maximum number of clusters
#' @return A vector of plot width and height
default_size <- function(max_cluster){
  unit_size <- 512*round(max_cluster/10)
  return(c(unit_size,unit_size))
}

#' @title get_layout
#' @description Get igraph layout by the number of nodes
#' @param subnet_condition A igraph sub-network
#' @import igraph
#' @return igraph network layout
get_layout <- function(subnet_condition) {
  if (length(V(subnet_condition)) >= 100){
    # layout_subnet <- layout_on_sphere(subnet_condition)
    # layout_subnet <- layout_with_kk(subnet_condition)
    # layout_subnet <- layout_with_lgl(subnet_condition)
    return(layout_with_fr(subnet_condition))
  } else {
    return(layout_nicely(subnet_condition))
  }
}

#' @title get_theme_blank
#' @description Predefined ggplot theme for removing ticks, titles and labels of X and Y axis
#' @import ggplot2
#' @return A ggplot theme
get_theme_blank <- function(){
  return(theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank()
  ))
}

#' @title length2
#' @description Customized function for vector length calculation
#' @param x A vector
#' @param na.rm If it is FALSE, no exclusion of NA values.
#' @return A vector length
length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)
}

#' @title get_stats
#' @description Calculate statistics of the given data for XINA network analysis
#' @param centrality_results Network centrality score data frame calculated by XINA network module
#' @param na.rm If it is FALSE, no exclusion of NA values.
#' @import plyr
#' @return A data frame containing statistics of XINA network centrality scores
get_stats <- function(centrality_results, na.rm=FALSE) {
  tgc <- ddply(data.frame(d=centrality_results$Score, r=centrality_results$Rank), c("r"), .drop=TRUE,
               .fun = function(xx, col) {
                 c(N    = length2(xx[[col]], na.rm=na.rm),
                   max  = max(xx[[col]], na.rm=na.rm),
                   min  = min(xx[[col]], na.rm=na.rm),
                   mean = mean   (xx[[col]], na.rm=na.rm)
                 )
               }, "d")
  return(tgc)
}

#' @title mutate_colors
#' @description 'mutate_colors' generates new color scheme for XINA clustering plot based on condition composition results (\link[XINA]{plot_condition_compositions}).
#' If any clusters have higher percentage than the 'threshold_percent', XINA will assign new colors in accordance to 'color_for_condition'.
#' If not, XINA will give 'gray' color or user-defined color via 'null_color' parameter.
#' @param condition_composition A data frame generated by \link[XINA]{plot_condition_compositions}
#' @param color_for_condition A vector like 'color_for_condition' of \link[XINA]{xina_clustering}
#' @param threshold_percent Default is 50.  The percentage threshold for giving new colors
#' @param null_color Default is 'gray'. This color is for clusters that are not biased to any of experimental conditions
#' @return A data frame containing statistics of XINA network centrality scores
#' @export
#' @examples
#' # load XINA example data
#' data(xina_example)
#'
#' # Plot condition composition pie-chart with default option
#' condition_composition <- plot_condition_compositions(example_clusters)
#' example_clusters$color_for_clusters <- mutate_colors(condition_composition,
#' example_clusters$color_for_condition)
#' plot_clusters(example_clusters, xval=c(0,2,6,12,24,48,72), xylab=FALSE)
#'
mutate_colors <- function(condition_composition, color_for_condition, null_color='gray', threshold_percent=50) {
  Cluster <- Condition <- NULL
  color_for_clusters <- c()
  conditions <- unique(condition_composition$Condition)
  cl=1
  for (cl in seq_len(max(condition_composition$Cluster))) {
    sub_df <- subset(condition_composition, Cluster==cl)
    color_for_cluster <- NULL
    for (condition in conditions) {
      sub_df2 <- subset(sub_df, Condition==condition)
      if (nrow(sub_df2) > 0)
        if (sub_df2$Percent_ratio >= threshold_percent)
          color_for_cluster <- as.character(color_for_condition[condition])
    }
    if (is.null(color_for_cluster)) {
      color_for_cluster <- null_color
    }
    color_for_clusters <- c(color_for_clusters, color_for_cluster)
  }
  return(color_for_clusters)
}

#' @title get_color_for_nodes
#' @description Pre-defined 30 colors
#' @return A vector for color code of XINA graphics
get_color_for_nodes <- function(){
  # http://colorbrewer2.org/#type=qualitative&scheme=Set3&n=12
  return(c('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854',
           '#ffd92f', '#e5c494', '#b3b3b3',
           '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
           '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a',
           '#ffff99', '#b15928',
           '#8dd3c7', '#ffffb3', '#bebada',
           '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5',
           '#d9d9d9', '#bc80bd', '#ccebc5', '#ffed6f'))
}

#' @title plot_NA
#' @description Draw NULL plot
#' @import graphics
#' @return a empty plot
plot_NA <- function() {
  plot(1, type="n", axes=FALSE, xlab="", ylab="", main="Not Available")
}

#' @title get_colors
#' @description Generate color series for XINA graphics
#' @param nClusters The number of clusters
#' @param set Pre-defined color series set
#' @param colorset manually defined color codes
#' @return A vector for color code of XINA graphics
get_colors <- function(nClusters, set='', colorset=NULL)
{
  if (is.null(colorset)) {
    if (set == '') {
      # Pastel tone
      colorset <- c("#FFB5E8", "#FF9CEE", "#FFCCF9", "#FCC2FF", "#F6A6FF",
                    "#B28DFF", "#C5A3FF", "#D5AAFF", "#ECD4FF", "#FBE4FF",
                    "#DCD3FF", "#A79AFF", "#B5B9FF", "#97A2FF", "#AFCBFF",
                    "#AFF8DB", "#C4FAF8", "#85E3FF", "#ACE7FF", "#6EB5FF",
                    "#BFFBCB", "#DBFFD6", "#BEBEBE", "#E7FFAC", "#FFFFD1",
                    "#FFC9DE", "#FFABAB", "#FFBEBC", "#FFCBC1", "#FFF5BA")
      # #FFFFD1 --> #d3d300
      N <- 30
    } else {
      # Original color
      colorset <- c("#8B008B", "#F1C40F", "#00008B","#1E90FF", "#104E8B",
                    "#FF0000", "#8B0000", "#9ACD32", "#008B00", "#00FFFF",
                    "#008B8B","#0000FF", "#8B0A50", "#FFFF00", "#8B7500",
                    "#00FF00", "#FF69B4", "#A020F0", "#FFA500", "#FF00FF")
      N <- 20
    }
  }
  if (nClusters>N){
    # color_for_graph <- colorRampPalette(colorset)(nClusters)
    color_for_graph <- c(rep(colorset, nClusters/length(colorset)),
                         colorset[seq_len(nClusters%%N)])
  } else {
    color_for_graph <- colorset
  }
  return(color_for_graph)
}

#' @title plot_condition_compositions
#' @description
#'  computes condition composition of the XINA clustering results and draws pie-charts.
#' @param clustering_result A list containing XINA clustering results.
#' See \link[XINA]{xina_clustering}
#' @param bullseye If it is TRUE, draw bullseye plot instead of the pie-chart. Default is FALSE
#' @param ggplot_theme This is ggplot theme to modify condition composition pie-chart and bulles eye plots.
#' @return A condition composition plot and a data frame containing condition compositions of the clusters
#' @import ggplot2
#' @import plyr
#' @import gridExtra
#' @export
#' @examples
#'
#' # load XINA example data
#' data(xina_example)
#'
#' # Plot condition composition pie-chart with default option
#' plot_condition_compositions(example_clusters)
#'
#' # Make a new color code for conditions
#' condition_colors <- c("tomato","steelblue1","gold")
#' names(condition_colors) <- example_clusters$condition
#' example_clusters$color_for_condition <- condition_colors
#'
#' # Draw condition composition pie-chart with the new color code
#' plot_condition_compositions(example_clusters)
#'
#' # Draw condition composition bullseye plot
#' plot_condition_compositions(example_clusters, bullseye = TRUE)
#'
plot_condition_compositions <- function(clustering_result, bullseye=FALSE, ggplot_theme=NULL)
{
  Cluster <- Condition <- Percent_ratio <- NULL
  # Collect cluster information by proteins
  clustering_dir <- clustering_result$out_dir
  nClusters <- clustering_result$nClusters
  max_cluster <- clustering_result$max_cluster
  data_column <- clustering_result$data_column
  column_numbers <- seq_len(length(data_column))
  super_ds <- clustering_result$clusters
  num_conditions <- length(unique(super_ds$Condition))
  # Color setting
  color_for_nodes <- clustering_result$color_for_condition
  if (length(color_for_nodes)==length(clustering_result$condition)) {
    names(color_for_nodes) <- clustering_result$condition
  } else {
    stop(paste("Length of color_for_nodes vector is not same to the number of conditions",
               length(color_for_nodes), length(clustering_result$condition)))
  }
  # Draw plots
  list_output <- list()
  out_summary <- data.frame()
  for (cn in seq_len(max_cluster))
  {
    cl_subset <- subset(super_ds, Cluster == cn)
    g_summary <- ddply(cl_subset, c("Condition"), .drop=TRUE,
                       .fun = function(xx, col) {
                         c(N    = length2(xx[[col]], na.rm=TRUE))},
                       c("Condition"))
    g_summary$Percent_ratio <- round((g_summary$N/sum(g_summary$N))*100, 2)
    out_summary <- rbind(out_summary,cbind(Cluster=cn, g_summary))
    plot_title <- paste("#",cn, " (n=",nrow(cl_subset),")", sep='')
    if (bullseye) {
      # Bullseye plot
      # http://ggplot2.tidyverse.org/reference/coord_polar.html
      list_output[[cn]] <- ggplot(cl_subset, aes(x = factor(1), fill = factor(Condition))) +
        geom_bar(width = 1) + coord_polar() +
        labs(title = plot_title) +
        scale_fill_manual(values=color_for_nodes) +
        get_theme_blank() +
        guides(fill=guide_legend(title=NULL))
    } else {
      # PieChart
      # http://www.sthda.com/english/wiki/ggplot2-pie-chart-quick-start-guide-r-software-and-data-visualization
      list_output[[cn]] <- ggplot(g_summary, aes(x="", y=Percent_ratio, fill=Condition)) +
        geom_bar(width = 1, stat = "identity") +
        coord_polar("y", start=0) +
        labs(title = plot_title) +
        scale_fill_manual(values=color_for_nodes) +
        get_theme_blank() +
        guides(fill=guide_legend(title=NULL))
    }
    if (!is.null(ggplot_theme)) {
      list_output[[cn]] <- list_output[[cn]] + ggplot_theme
    }
  }
  do.call(grid.arrange,list_output)
  return(out_summary)
}

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

#' @title plot_clusters
#' @description Draw all the clustering results.
#' 'plot_clusters' draws two plots, scaled and unscaled line graphs.
#' Scaled graphs have same y limits that are 0 to 1 by default,
#' but can be changed via 'y_lim' parameter.
#' @param clustering_result A list containing XINA clustering results.
#' See \link[XINA]{xina_clustering}
#' @param y_lim Y axis limit. If you set y_lim=c(0,1),
#' 'plot_clusters' will plot line graphs scaled from 0 to 1 in y-axis
#' Default is NULL, which means unscaled line graphs.
#' @param xval Change X axis values and labels.
#' Default is data_column of the clustering result list
#' @param xylab If it is FALSE, x and y labels will be blank.
#' If it is TRUE (defualt), x and y labels will be shown.
#' @param ggplot_theme This is ggplot theme to modify XINA clustering plot.
#' @return Line graphs of all the clusters
#' @import ggplot2
#' @export
#' @examples
#' library(ggplot2)
#'
#' # load XINA example data
#' data(xina_example)
#'
#' # Draw clustering plots
#' plot_clusters(example_clusters)
#'
#' # Apply theme to the clustering plot
#' theme1 <- theme(title=element_text(size=8, face='bold'),
#' axis.text.x = element_text(size=7),
#' axis.text.y = element_blank(),
#' axis.ticks.x = element_blank(),
#' axis.ticks.y = element_blank(),
#' axis.title.x = element_blank(),
#' axis.title.y = element_blank())
#' plot_clusters(example_clusters, ggplot_theme=theme1)
#'
plot_clusters <- function(clustering_result, y_lim=NULL, xval=NULL, xylab=TRUE, ggplot_theme=NULL)
{
  Cluster <- indviduals <- x <- y <- NULL
  clustering_dir <- clustering_result$out_dir
  nClusters <- clustering_result$nClusters
  max_cluster <- clustering_result$max_cluster
  data_column <- clustering_result$data_column
  data <- clustering_result$clusters
  # X axis set up
  if (is.null(xval)) {
    column_numbers <- seq_len(length(data_column))
    xtick_label <- data_column
  } else {
    column_numbers <- xval
    xtick_label <- xval
  }
  # Plot parameter setting
  color_for_graph <- clustering_result$color_for_clusters
  if (length(color_for_graph)<max_cluster) {
    stop(paste("Length of color_for_clusters vector should be bigger than max_cluster",
               length(clustering_result$color_for_clusters), max_cluster))
  }
  # Draw line graphs and save it as a PNG image
  scaled_plots <- list()
  plot_out <- list()
  for (i in seq_len(max_cluster)) {
    sub_data<-as.matrix(subset(data, Cluster %in% i, select = data_column))
    nlines <- nrow(sub_data)
    if (nlines > 0) {
      data_Sds<-cbind(data.frame(expand.grid(x = column_numbers, indviduals = seq_len(nlines))),
                      y=as.numeric(t(sub_data[, data_column])))
      num_prots <- nrow(data_Sds)/length(data_column)
      plt_title <- paste("#",i, " (n=",num_prots,")", sep='')
      plot_out[[i]]<-ggplot(data_Sds, aes(x=x, y=y, group=indviduals)) +
        geom_line(colour = color_for_graph[i]) +
        geom_hline(yintercept=1/length(data_column), color="red", linetype="dashed") +
        labs(x="", y="", title=plt_title) +
        scale_x_continuous(breaks=column_numbers, labels=xtick_label)
      if(!xylab){
        plot_out[[i]] <- plot_out[[i]] + get_theme_blank()
      }
      if (!is.null(y_lim)) {
        plot_out[[i]] <- plot_out[[i]] + ylim(y_lim)
      }
      if (!is.null(ggplot_theme)) {
        plot_out[[i]] <- plot_out[[i]] + ggplot_theme
      }
    } else {
      plot_NA()
    }
  }
  do.call(gridExtra::grid.arrange,plot_out)
  # return(plot_out)
}

#' @title plot_clusters_all
#' @description Draw line graphs of all the proteins in the given dataset
#' @param clustering_result A list containing XINA clustering results.
#' See \link[XINA]{xina_clustering}
#' @param selected_condition A condition name to draw the kinetics plot
#' @return a list containing clustering results and pdf file
#' containing a BIC plot in current working directory.
#' @import ggplot2
#' @export
#' @examples
#'
#' # load XINA example data
#' data(xina_example)
#'
#' # Plot kinetics of all the proteins in Control
#' plot_clusters_all(example_clusters, selected_condition="Control")
#'
#' # Plot kinetics of all the proteins in Stimulus1
#' plot_clusters_all(example_clusters, selected_condition="Stimulus1")
#'
#' # Plot kinetics of all the proteins in Stimulus2
#' plot_clusters_all(example_clusters, selected_condition="Stimulus2")
#'
#' # Plot kinetics of all the proteins in three data
#' plot_clusters_all(example_clusters)
#'
plot_clusters_all <- function(clustering_result, selected_condition=NULL)
{
  condition <- x <- y <- z_levels_Sds <- NULL
  data_column <- clustering_result$data_column
  data <- clustering_result$clusters
  column_numbers <- seq_len(length(data_column))
  color_for_condition <- clustering_result$color_for_condition
  all_prots <- data.frame(expand.grid(x = seq_len(length(data_column)), y=seq_len(nrow(data))),
                          condition = as.vector(data[,"Condition"]),
                          z_levels_Sds=as.numeric(t(data[data_column])))
  if (!is.null(selected_condition)) {
    all_prots <- subset(all_prots, condition==selected_condition)
    title2show <- paste("Unclustered profiles of",selected_condition)
  } else {
    title2show <- "Unclustered profiles of all"
  }
  ggplot(all_prots, aes(x=x, y=z_levels_Sds, condition=y, colour=condition)) +
    geom_line() +
    scale_color_manual(values=color_for_condition) +
    geom_hline(yintercept=1/max(column_numbers), color="red", linetype="dashed")+
    labs(x = "", y = paste("Normalized intensity by",clustering_result$norm_method)) +
    labs(title =  title2show) +
    scale_x_continuous(breaks=column_numbers, labels=data_column)
}


#' @title get_condition_biased_comigrations
#' @description get comigrations that at least one biased cluster is involved in.
#' Biased clusters are defined by
#' @param clustering_result A list containing XINA clustering results.
#' See \link[XINA]{xina_clustering}
#' @param count_table A data frame generated by using \link[plyr]{count}.
#' If count_table is NULL (by default), XINA will consider all the comigrations.
#' @param selected_conditions A vector of condition names used in XINA clustering results.
#' The number of selected conditions should be at least two.
#' @param condition_composition The resulting data frame of 'plot_condition_compositions'.
#' See \link[XINA]{plot_condition_compositions}.
#' @param threshold_percent Default is 50.  The percentage threshold for finding condition-biased clusters
#' @param color_for_null A color for non-condition-biased comigrations. Default is 'gray'
#' @param color_for_highly_matched A color for comigrations that are involved with more than two condition-biased clusters. Default is 'red4'
#' @param cex Size of cluster number on block axis. Default if 0.7. See \link[alluvial]{alluvial}.
#' @param alpha Transparency of alluvia colors. Default is 0.3. See \link[alluvial]{alluvial}.
#' @return An alluvial plot displaying comigrations and the data frame containing condition-biased comigrations.
#' @export
#' @examples
#'
#' # load XINA example data
#' data(xina_example)
#'
#' # get a vector of experimental conditions analyzed in the clustering results
#' conditions <- as.vector(example_clusters$condition)
#'
#' # get condition composition information
#' condition_composition <- plot_condition_compositions(example_clusters)
#'
#' comigrations_size10 <- alluvial_enriched(example_clusters, conditions, comigration_size=10)
#' # Finding condition-biased comigrations by 50% threshold
#' condition_biased_comigrations <-
#' get_condition_biased_comigrations(clustering_result=example_clusters,
#' count_table=comigrations_size10, selected_conditions=conditions,
#' condition_composition=condition_composition)
#'
#' # Finding condition-biased comigrations by 70% threshold
#' condition_biased_comigrations <-
#' get_condition_biased_comigrations(clustering_result=example_clusters,
#' count_table=comigrations_size10, selected_conditions=conditions,
#' condition_composition=condition_composition,
#' threshold_percent=70)
#'
get_condition_biased_comigrations <- function(clustering_result, count_table=NULL,
                                              selected_conditions, condition_composition, threshold_percent=50,
                                              color_for_null='gray', color_for_highly_matched='red4',
                                              cex=0.7, alpha=0.3) {
  Condition <- Percent_ratio <- NULL
  if (is.null(count_table)) {
    count_table <- generate_count_table(clustering_result, selected_conditions, 0)
  }
  ratio_threshold <- threshold_percent
  biased_clusters <- list()
  for (condition in selected_conditions) {
    tmp <- subset(condition_composition, Percent_ratio>=ratio_threshold & Condition==condition)
    biased_clusters[[condition]] <- tmp$Cluster
  }
  alluvia_colors <- c()
  for (i in seq_len(nrow(count_table))) {
    matched_conditions <- c()
    for (condition in selected_conditions) {
      cl <- count_table[i,condition]
      if (!is.na(match(cl, biased_clusters[[condition]]))) {
        matched_conditions <- c(matched_conditions, condition)
      }
    }
    if (length(matched_conditions)==1) {
      alluvia_colors <- c(alluvia_colors, clustering_result$color_for_condition[matched_conditions])
    } else if (length(matched_conditions)>1) {
      alluvia_colors <- c(alluvia_colors, color_for_highly_matched)
    } else {
      alluvia_colors <- c(alluvia_colors, color_for_null)
    }
  }
  cnt_table_4draw <- count_table[!is.na(match(alluvia_colors, c(clustering_result$color_for_condition, color_for_highly_matched))),]
  alluvia_colors <- alluvia_colors[!is.na(match(alluvia_colors, c(clustering_result$color_for_condition, color_for_highly_matched)))]
  if (nrow(cnt_table_4draw)>0){
    return(draw_alluvial_plot(clustering_result, selected_conditions, cnt_table_4draw,
                              alluvia_colors = alluvia_colors, cex=cex, alpha=alpha))
  } else {
    print("No biased comigration was found")
  }
}

#' @title get_comigrations_by_name
#' @description 'get_comigrations_by_name' finds proteins comigrated with the given proteins
#' @param clustering_result A list containing XINA clustering results. See \link[XINA]{xina_clustering}
#' @param selected_conditions A vector of condition names used in XINA clustering results.
#' The number of selected conditions should be at least two.
#' @param protein_list A vector containing gene names.
#' @param cex Size of cluster number on block axis. Default if 0.7. See \link[alluvial]{alluvial}
#' @param alpha Transparency of alluvia colors. Default is 0.3. See \link[alluvial]{alluvial}
#' @return An alluvial plot displaying comigrations and the data frame containing comigrations of the input proteins
#' @export
#' @examples
#'
#' # load XINA example data
#' data(xina_example)
#'
#' # the clustering result table
#' all_proteins  <- as.character(example_clusters$aligned$`Gene name`)
#' # get a vector of experimental conditions analyzed in the clustering results
#' classes <- as.vector(example_clusters$condition)
#'
#' comigrated_prots_all <- get_comigrations_by_name(example_clusters, classes, all_proteins[1:3])
#'
get_comigrations_by_name <- function(clustering_result, selected_conditions, protein_list, cex=0.7, alpha=0.3){
  matched_comigrations_all <- data.frame()
  comigrated_prots_all <- list()
  count_table <- generate_count_table(clustering_result, selected_conditions, 0)
  found_proteins <- c()
  for (prot_name in protein_list) {
    cl_matched <- clustering_result$aligned[clustering_result$aligned$"Gene name"==prot_name, selected_conditions]
    comigrations_matched <- count_table
    comigrated_prots <- clustering_result$aligned
    flag_found <- TRUE
    if (length(cl_matched)==0) {
      flag_found <- FALSE
    } else if (is.na(sum(as.numeric(cl_matched)))) {
      flag_found <- FALSE
    } else if (sum(as.numeric(cl_matched))==0) {
      flag_found <- FALSE
    }
    if (flag_found){
      if (length(selected_conditions)==1){
        names(cl_matched) <- selected_conditions
      }
      for (condition in selected_conditions){
        comigrations_matched <- subset(comigrations_matched, comigrations_matched[condition]==as.integer(cl_matched[condition]))
        comigrated_prots <- subset(comigrated_prots, comigrated_prots[condition]==as.integer(cl_matched[condition]))
      }
      found_proteins <- c(found_proteins, prot_name)
      comigrated_prots_all[[prot_name]] <- comigrated_prots
      if (nrow(matched_comigrations_all)==0) {
        matched_comigrations_all <- comigrations_matched
      } else {
        matched_comigrations_all <- rbind(matched_comigrations_all,comigrations_matched)
      }
    } else {
      print(paste(prot_name," is not found in the your 'selected_conditions':",selected_conditions))
      comigrated_prots_all[[prot_name]] <- NULL
    }
  }
  rownames(matched_comigrations_all) <- found_proteins
  if (length(selected_conditions)>1){
    comigrated_prots_all[["comigrations"]] <- draw_alluvial_plot(clustering_result, selected_conditions, matched_comigrations_all,
                                                                 cex=cex, alpha=alpha)
  } else {
    comigrated_prots_all[["comigrations"]] <- matched_comigrations_all
  }
  return(comigrated_prots_all)
}

#' @title draw_alluvial_plot
#' @description 'draw_alluvial_plot' draw a alluvial plot
#' @param clustering_result A list containing XINA clustering results.
#' See \link[XINA]{xina_clustering}.
#' @param selected_conditions A vector of condition names used in XINA clustering results.
#' The number of selected conditions should be at least two.
#' @param count_table A data frame generated by using \link[plyr]{count}.
#' @param alluvia_colors A vector containing the user-defined colors for each alluvium.
#' @param cex Size of cluster number on block axis. Default if 0.7. See \link[alluvial]{alluvial}.
#' @param alpha Transparency of alluvia colors. Default is 0.3. See \link[alluvial]{alluvial}.
#' @return An alluvial plot displaying comigrations and the data frame containing the input count_table with colors.
#' @import alluvial
#' @export
#' @examples
#'
#' # load XINA example data
#' data(xina_example)
#'
#' # get a vector of experimental conditions analyzed in the clustering results
#' classes <- as.vector(example_clusters$condition)
#'
#' comigrations_size_over5 <- alluvial_enriched(example_clusters, classes, comigration_size=5)
#' draw_alluvial_plot(example_clusters, classes, comigrations_size_over5)
#'
draw_alluvial_plot <- function(clustering_result, selected_conditions, count_table,
                               alluvia_colors=NULL, cex=0.7, alpha=0.3){
  # Set color codes
  if (is.null(alluvia_colors)) {
    alluvia_colors <- c()
    color_codes <- get_colors(clustering_result$nClusters, set='vivid')
    for (i in seq_len(nrow(count_table))){
      cn <- as.numeric(as.vector(count_table[i,1]))
      if (cn == 0){
        alluvia_colors <- c(alluvia_colors, "#BEBEBE")
      } else {
        alluvia_colors <- c(alluvia_colors, color_codes[cn])
      }
    }
  } else {
    if (length(alluvia_colors) > nrow(count_table)){
      alluvia_colors <- alluvia_colors[seq_len(nrow(count_table))]
    } else if (length(alluvia_colors) < nrow(count_table)){
      difference <- nrow(count_table) - length(alluvia_colors)
      for (i in seq_len(difference)) {
        alluvia_colors <- c(alluvia_colors, 'gray')
      }
      warning(paste(difference, "comigrations were colored by gray because alluvial_colors is not enough to color all the comigrations"))
    }
  }
  # For the dimension error when there is only one row in count_table
  if (nrow(count_table)==0){
    stop("No biased comigration was found")
  } else if (nrow(count_table)==1){
    Comigration_size <- c(count_table$Comigration_size,0)
    cols <- c(alluvia_colors,alluvia_colors)
    bds <- c(alluvia_colors,alluvia_colors)
  } else {
    Comigration_size <- count_table$Comigration_size
    cols <- alluvia_colors
    bds <- alluvia_colors
  }
  # Draw an alluvial plot
  # count_table[,1:length(selected_conditions)]
  alluvial(count_table[,selected_conditions],
           freq=Comigration_size,
           col=cols,
           border=bds,
           cex=cex,
           alpha=alpha)
  if (is.null(count_table$"Alluvia_color")) {
    return(cbind(count_table, Alluvia_color=alluvia_colors))
  } else {
    count_table$"Alluvia_color" <- alluvia_colors
    return(count_table)
  }
}


#' @title alluvial_enrichment_tests
#' @description Fisher's exact test to calculate the significance over all comigrations.
#' The following 2x2 table was used to calculate p-value from Fisher's exact test.
#' To evaluate significance of comigrated proteins from cluster #1 in control to cluster #2 in test condition,
#'      \tabular{rll}{
#'            \tab \strong{cluster #1 in control} \tab \strong{other clusters in control}\cr
#'       \strong{cluster #2 in test}      \tab 65 (TP) \tab 175 (FP)\cr
#'       \strong{other clusters in test}  \tab 35 (FN) \tab 979 (TN)\cr
#'      }
#' 'alluvial_enrichment_tests' also provides another statistical methods including Hypergeometric test and Chi-square test.
#' @param count_table A data frame generated by using \link[plyr]{count}.
#' @param c1 A selected cluster in the first condition.
#' @param c2 A selected cluster in the second condition.
#' @param non_cluster The cluster number for proteins that were not detected in a specific sample. Default is 0.
#' @param test_type Enrichment test type.
#' 'fisher' = Fisher's exact test, 'hyper' = Hypergeometric test, 'chisq' = Chi-square test
#' @return P-value of comigration enrichment test and 2x2 table information
alluvial_enrichment_tests <- function(count_table, c1, c2, non_cluster=0, test_type='fisher'){
  matched_in_list <- sum(subset(count_table, count_table[,1]==c1 & count_table[,2]==c2)$Comigration_size)
  non_matched_in_list <- sum(subset(count_table, count_table[,1]==c1 & count_table[,2]!=c2)$Comigration_size)
  matched_in_all <- sum(subset(count_table, count_table[,1]!=c1 & count_table[,2]==c2)$Comigration_size)
  non_matched_in_all <- sum(subset(count_table, count_table[,1]!=c1 & count_table[,2]!=c2)$Comigration_size)
  data_table <- c(matched_in_list, non_matched_in_list, matched_in_all, non_matched_in_all)
  counts = (matrix(data = data_table, nrow = 2))
  pval <- 1
  if (test_type=='fisher') {
    pval <- stats::fisher.test(counts)$p.value
  } else if (test_type=='hyper') {
    pval <- stats::dhyper(matched_in_list, matched_in_all, non_matched_in_all, matched_in_list+non_matched_in_list)
  } else if (test_type=='chisq') {
    pval <-  stats::chisq.test(counts)$p.value
  } else {
    stop("You chose wrong test_type that should be one of c('fisher', 'hyper', 'chisq')")
  }
  return_table <- c(data_table, pval)
  names(return_table) <- c("TP","FP","FN","TN","Pvalue")
  return(return_table)
}

#' @title generate_count_table
#' @description Count the number of comigrated proteins using \link[plyr]{count}
#' @param clustering_result A list containing XINA clustering results. See \link[XINA]{xina_clustering}
#' @param selected_conditions A vector of condition names used in XINA clustering results.
#' @param comigration_size The number of proteins comigrated together in the selected conditions of XINA clustering results. Default is 0.
#' @return A data frame containing comigrations.
generate_count_table <- function(clustering_result, selected_conditions, comigration_size) {
  # Generate count table
  aligned_clusters <- clustering_result$aligned
  selected_data <- data.frame(Gene_name=aligned_clusters$`Gene name`)
  for (i in selected_conditions){
    selected_data <- cbind(selected_data, as.integer(aligned_clusters[i][[1]]))
  }
  colnames(selected_data) <- c("Gene_name",selected_conditions)
  count_table <- plyr::count(selected_data[selected_conditions])
  count_table$"RowNum" <- seq_len(nrow(count_table))
  na_cluster <- count_table
  for (i in selected_conditions){
    if (nrow(na_cluster)>0) {
      na_cluster <- subset(na_cluster, na_cluster[i]==0)
    }
  }
  return_table <- count_table[!(count_table$RowNum %in% na_cluster$RowNum),]
  colnames(return_table) <- c(selected_conditions, "Comigration_size","RowNum")
  return(return_table)
}

#' @title alluvial_enriched
#' @description 'alluvial_enriched’ draws an alluvial plot and finds comigrated proteins.
#' The comigration is a group of proteins that show the same expression pattern,
#' classified and evaluated by XINA clustering, in at least two conditions.
#' XINA can reduce the dataset complexity by filtering based on
#' the number of comigrated proteins (size, ’comigration_size’ parameter) and
#' perform an enrichment test (P-value of Fisher’s exact test, ’pval_threshold’)
#' to determine significance of enriched comigrations.
#' The Fisher’s exact test can only be done for two conditions at a time.
#' The following 2x2 table was used to calculate the P-value from the Fisher’s exact test.
#' To evaluate significance of co-migrated proteins from cluster #1 in control to cluster #2 in test group,
#'      \tabular{rll}{
#'       - \tab cluster #1 in control \tab other clusters in control\cr
#'       cluster #2 in test      \tab 65 (TP) \tab 175 (FP)\cr
#'       other clusters in test  \tab 35 (FN) \tab 979 (TN)\cr
#'      }
#' @param clustering_result A list containing XINA clustering results. See \link[XINA]{xina_clustering}
#' @param selected_conditions A vector of condition names used in XINA clustering results.
#' The number of selected conditions should be at least two.
#' @param comigration_size The number of proteins comigrated together in the selected conditions of XINA clustering results.
#' Default is 0
#' @param pval_threshold This option is avaiable only when you selected two conditions for comigration search.
#' @param pval_method Method for p-value adjustment. See \link[stats]{p.adjust}
#' @param cex Scaling of fonts of category labels. Default if 0.7. See \link[alluvial]{alluvial}
#' @param alpha Transparency of the stripes. Default if 0.3. See \link[alluvial]{alluvial}
#' @return A data frame containing comigrations and an alluvial plot showing comigrations
#' @export
#' @examples
#'
#' # load XINA example data
#' data(xina_example)
#'
#' # Get the experimental conditions in the example data
#' classes <- as.vector(example_clusters$condition)
#'
#' # Get comigrations without any thresholds
#' all_comigrations <- alluvial_enriched(example_clusters, classes)
#'
#' # Get comigrations that have >= 5 size (the number of comigrated proteins)
#' all_cor_enriched <- alluvial_enriched(example_clusters, classes, comigration_size=5)
#'
#' # Get all the comigrations between Control and Stimulus1
#' comigrations_Control_Stimulus1 <- alluvial_enriched(example_clusters,
#' c(classes[1],classes[2]))
#'
#' # Get comigrations between Control and Stimulus1, that have >=5 size
#' comigrations_Control_Stimulus1_over5 <- alluvial_enriched(example_clusters,
#' c(classes[1],classes[2]), comigration_size=5)
#'
#' # Get comigrations between Control and Stimulus1,
#' # that have >= 5 size and enrichment FDR <= 0.01
#' comigrations_Control_Stimulus1_pval0.01_size5 <- alluvial_enriched(example_clusters,
#' c(classes[1],classes[2]), comigration_size=5, pval_threshold=0.01)
#'
#' # Get  comigrations between Control and Stimulus1,
#' # that have >= 5 size and enrichment Benjamini & Yekutieli <= 0.01
#' comigrations_Control_Stimulus1_BY0.01_size5 <- alluvial_enriched(example_clusters,
#' c(classes[1],classes[2]), comigration_size=5, pval_threshold=0.01, pval_method="BY")
#'
alluvial_enriched <- function(clustering_result, selected_conditions, comigration_size=0,
                              pval_threshold=1, pval_method='fdr', cex=0.7, alpha=0.3) {
  Comigration_size <- Pvalue.adjusted <- NULL
  # For debugging
  # selected_conditions <- c(classes[1],classes[2],classes[3])
  # https://cran.r-project.org/web/packages/alluvial/vignettes/alluvial.html#changing-layers
  if (length(selected_conditions)<2){
    stop("'selected_conditions' requires at least two")
  } else if (length(selected_conditions)>2) {
    warning("length(selected_conditions) > 2, so XINA can't apply the enrichment filter
            Can't apply the enrichment filter, so pval_threshold is ignored")
    pval_threshold <- 1
  }
  # generate a count_table
  count_table <- generate_count_table(clustering_result, selected_conditions, 0)

  # Do enrichment test
  tp <- c()
  fp <- c()
  fn <- c()
  tn <- c()
  p.values <- c()
  if (length(selected_conditions)==2){
    for (j in seq_len(nrow(count_table))) {
      c1 <- count_table[j,1]
      c2 <- count_table[j,2]
      test_result <- alluvial_enrichment_tests(count_table,c1,c2)
      tp <- c(tp, test_result["TP"])
      fp <- c(fp, test_result["FP"])
      fn <- c(fn, test_result["FN"])
      tn <- c(tn, test_result["TN"])
      p.values <- c(p.values, test_result["Pvalue"])
    }
    pval_adjusted <- stats::p.adjust(p.values,method=pval_method)
  } else {
    p.values <- rep(NA,nrow(count_table))
    pval_adjusted <- rep(NA,nrow(count_table))
    tp <- rep(NA,nrow(count_table))
    fp <- rep(NA,nrow(count_table))
    fn <- rep(NA,nrow(count_table))
    tn <- rep(NA,nrow(count_table))
  }
  count_table_all <- cbind(count_table,
                           PValue=p.values,
                           Pvalue.adjusted=pval_adjusted,
                           TP=tp,FP=fp,FN=fn,TN=tn)
  if (pval_threshold == 1) {
    count_table <- subset(count_table_all, Comigration_size>=comigration_size)
  } else {
    count_table <- subset(count_table_all, Pvalue.adjusted<=pval_threshold & Comigration_size>=comigration_size)
  }
  if (nrow(count_table)>0){
    count_table_colored <- draw_alluvial_plot(clustering_result, selected_conditions, count_table, cex=cex, alpha=alpha)
  } else {
    print("Can't draw alluvial plot because none of features passed the thresholds")
  }
  return(count_table_colored)
}


#' @title calculate_centrality_scores
#' @description 'calculate_centrality_scores' computes network centrality scores
#' @param net protein-protein interaction network of igraph
#' @param centrality_type the maximum number of clusters
#' @import igraph
#' @return A vector of network centrality scores
calculate_centrality_scores <- function(net, centrality_type="Degree") {
  if (centrality_type == "Degree") {
    centrality_score <- degree(net, mode="all")
  } else if (centrality_type == "Eigenvector") {
    centrality_score <- eigen_centrality(net, directed=FALSE, weights=NA)$vector
  } else if (centrality_type == "Hub") {
    centrality_score <- hub_score(net, weights=NA)$vector
  } else if (centrality_type == "Authority") {
    centrality_score <- authority_score(net, weights=NA)$vector
  } else if (centrality_type == "Closeness") {
    centrality_score <- closeness(net, mode="all", weights=NA)
  } else if (centrality_type == "Betweenness") {
    centrality_score <- betweenness(net, directed=FALSE, weights=NA)
  } else {
    stop("'centrality_type' should be one of c('Degree', 'Eigenvector', 'Hub', 'Authority',
         'Closeness', 'Betweenness')")
  }
  return(as.numeric(as.character(centrality_score)))
  }

#' @title rank_centrality
#' @description Give ranks based on network centrality scores
#' @param centrality_score Network centrality score matrix
#' @param type Network centrality score type, such as 'Eigenvector'
#' @param num_breaks The number of ranks
#' @return A vector containing ranks
rank_centrality <- function(centrality_score, type, num_breaks=5){
  if (type == "Closeness") {
    rank <- as.numeric(cut(rank(centrality_score), breaks=num_breaks))
  } else {
    if (length(centrality_score)<num_breaks){
      print("Since the number of centrality scores are less than 'num_breaks', so rank is not calculated")
      return(rep(1,length(centrality_score)))
    }
    rank <- as.numeric(cut(centrality_score, breaks=num_breaks))
    if (length(unique(rank))==1){
      rank <- rep(1,length(rank))
    }
  }
  return(rank)
}

#' @title xina_enrichment
#' @description xina_enrichment conducts functional enrichment tests using gene ontology or KEGG pathway terms for a given protein list
#' @param string_db STRINGdb object
#' @param protein_list A vector of gene names to draw protein-protein interaction network.
#' @param enrichment_type A functional annotation for the enrichment test.
#' 'enrichment_type' should be one of 'GO' and 'KEGG',
#' @param pval_threshold P-value threshold to get significantly enriched terms from the given proteins
#' @param methodMT Method for p-value adjustment. See \link[STRINGdb]{get_enrichment}. Default is 'fdr'.
#' @return A list of data frames containing enrichment results
#' @export
#' @import STRINGdb
#' @examples
#' \dontrun{
#' library(STRINGdb)
#' library(Biobase)
#'
#' # load XINA example data
#' data(xina_example)
#'
#' # Get STRING database for protein-protein intereaction information
#' string_db <- STRINGdb$new( version="10", species=9606, score_threshold=0, input_directory="" )
#' string_db
#'
#' # XINA analysis with STRING DB
#' xina_result <- xina_analysis(example_clusters, string_db)
#'
#' # Select proteins that showed cluster #1 in the Stimulus2 condition
#' subgroup <- subset(example_clusters$aligned, Stimulus2==1)
#' protein_list <- as.vector(subgroup$`Gene name`)
#'
#' # Enrichment test using KEGG pathway terms that have adjuseted p-value less than 0.1
#' kegg_enriched <- xina_enrichment(string_db, protein_list,
#' enrichment_type = "KEGG", pval_threshold=0.1)
#' plot_enrichment_results(kegg_enriched$KEGG, num_terms=10)
#'
#' # Enrichment test using GO terms that have adjuseted p-value less than 0.1
#' go_enriched <- xina_enrichment(string_db, protein_list,
#' enrichment_type = "GO", pval_threshold=0.1)
#' plot_enrichment_results(go_enriched$Component, num_terms=10)
#' }
#'
xina_enrichment <- function(string_db, protein_list, enrichment_type="GO",
                            pval_threshold=0.05, methodMT='fdr') {
  STRING_id <- pvalue_fdr <- NULL
  string_id <- string_db$map(data.frame(Accession=protein_list), "Accession")
  string_id <- subset(string_id, !is.na(STRING_id))

  if (enrichment_type=="GO"){
    categories <- c("Process","Function","Component")
  } else if (enrichment_type=="KEGG"){
    categories <- c(enrichment_type)
  } else {
    stop("'enrichment_type' should be one of c(GO, KEGG)")
  }
  return_list <- list()
  for (cat in categories){
    enrich_result <- string_db$get_enrichment(string_id$STRING_id, category=cat, methodMT=methodMT, iea=TRUE)
    return_list[[cat]] <- subset(enrich_result, pvalue_fdr<pval_threshold)
  }
  return(return_list)
}

#' @title xina_analysis
#' @description xina_analysis is to analyze protein-protein interaction(PPI) networks using STRINGdb and igraph R package. This module computes PPI networks within each XINA clusters.
#' @param clustering_result A list containing XINA clustering results. See \link[XINA]{xina_clustering}
#' @param ppi_db STRINGdb object
#' @param is_stringdb If it is TRUE (default), XINA will process 'ppi_db' as STRINGdb,
#' but it is FALSE, XINA will accepts your 'ppi_db' as it is.
#' You can make your own igraph network using customized PPI information instead of STRINGdb.
#' @param flag_simplify If it is TRUE (default), XINA will exclude unconnected proteins
#' @param node_shape You can choose node shape. Default is "sphere".  See \link[igraph]{shapes}
#' @param num_clusters_in_row The number of clusters in a row on the XINA network plot. Default is 5.
#' @param img_size Set the image size. For width=1000 and height=1500, it is img_size=c(1000,1500).
#' @param img_qual Set the image resolution. Default is 300.
#' @return A PNG file (XINA_Cluster_Networks.png) displaying PPI network plots of all the clusters
#' and a list containing XINA network analysis results.
#'      \tabular{rl}{
#'       \strong{Item} \tab \strong{Description}\cr
#'       All_network \tab PPI network of all the input proteins\cr
#'       Sub_network \tab A list containing PPI networks of each clusters\cr
#'       Data \tab XINA clustering results. See \link[XINA]{xina_clustering}\cr
#'       Nodes \tab A list of proteins in each cluster\cr
#'       Conditions \tab A list of experimental condition of proteins in each cluster\cr
#'       Titles \tab  A list of plot titles for XINA plotting\cr
#'       out_dir \tab A directory path storing XINA network analysis results\cr
#'       is_stringdb \tab False = different PPI DB and TRUE = STRING DB\cr
#'      }
#' @import igraph
#' @import grDevices
#' @import graphics
#' @export
#' @examples
#' \dontrun{
#' # load XINA example data
#' data(xina_example)
#'
#' # use the following code for utilizing up-to-date STRING DB
#' tax_id <- 9606  # for human
#' # tax_id <- 10090  # for mouse
#' library(STRINGdb)
#' library(igraph)
#' string_db <- STRINGdb$new( version='10', species=tax_id, score_threshold=0, input_directory='' )
#' string_db
#' xina_result <- xina_analysis(example_clusters, string_db, flag_simplify=FALSE)
#'
#' # Run XINA with a protein-protein interaction edgelist
#' data(HPRD)
#' net_all <- simplify(graph_from_data_frame(d=hprd_ppi, directed=FALSE),
#' remove.multiple = FALSE, remove.loops = TRUE)
#' xina_result <- xina_analysis(example_clusters, net_all, is_stringdb=FALSE, flag_simplify=FALSE)
#' }
#'
xina_analysis <- function(clustering_result, ppi_db, is_stringdb=TRUE, flag_simplify=TRUE,
                          node_shape="sphere", num_clusters_in_row=5, img_size=NULL, img_qual=300) {
  STRING_id <- NULL
  # Collect cluster information by proteins
  na_dir <- clustering_result$out_dir
  nClusters <- clustering_result$nClusters
  data_column <- clustering_result$data_column
  column_numbers <- seq_len(length(data_column))
  super_ds <- clustering_result$clusters
  num_conditions <- length(unique(super_ds$Condition))
  uniq_condition <- unique(super_ds$Condition)
  max_cluster <- clustering_result$max_cluster
  num_clusters_in_col <- as.integer(max_cluster/num_clusters_in_row) + if(max_cluster%%num_clusters_in_row>0){1}else{0}
  if (!is.vector(img_size)){
    unit_size <- 750
    img_size <- c(unit_size*num_clusters_in_row,unit_size*num_clusters_in_col)
  }
  # get required colors
  color_for_graph <- clustering_result$color_for_clusters
  color_for_nodes <- clustering_result$color_for_condition
  # analyze STRINGdb
  if (is_stringdb){
    superset_string <- ppi_db$map(super_ds, "Accession")
    colnames(superset_string) <- c("Accession", data_column,
                                   "Description", "Condition",
                                   "Key", "Cluster", "STRING_id")
    net_stringDB <- ppi_db$graph
    vertices <- as.vector(V(net_stringDB)$name)
    # Get protein lists that are matched to the PPI network
    prots <- unique(superset_string$STRING_id)
    nodes <- c()
    for (k in seq_len(length(prots))) {
      nodes <- c(nodes, match(prots[k], vertices))
    }
    # Get a subnetwork using the matched proteins
    subnet <- induced_subgraph(net_stringDB, nodes[!is.na(nodes)])
    SIDs <- V(subnet)$name
    sub_acc <- c()
    sub_desc <- c()
    for (k in seq_len(length(SIDs))) {
      sub_tmp <- subset(superset_string, STRING_id==SIDs[k])
      sub_acc <- c(sub_acc, as.vector(sub_tmp$Accession)[1])
      sub_desc <- c(sub_desc, as.vector(sub_tmp$Description)[1])
    }
    V(subnet)$name <- sub_acc
    V(subnet)$desc <- sub_desc
    net_all <- subnet
    # analyze user-defined PPI database
  } else {
    superset_string <- super_ds
    net_all <- ppi_db
  }
  vertices <- as.vector(V(net_all)$name)
  # Draw all the XINA plots into one page
  out <- paste(na_dir,"/","XINA_Cluster_Networks.png",sep="")
  png(filename=out, width=img_size[1], height=img_size[2], res=img_qual)
  par(mfrow = c(num_clusters_in_col, num_clusters_in_row))
  par(mar=c(1,1,1,1))
  list_nodes <- list()
  list_conditions <- list()
  list_subnets <- list()
  titles <- c()
  # Draw networks
  for (i in seq_len(max_cluster)) {
    prots_clustered <- subset(super_ds, super_ds$Cluster==as.character(i))
    if (nrow(prots_clustered)>0) {
      prots <- toupper(prots_clustered$Accession)
      conditions <- prots_clustered$Condition
      # Get protein lists that are matched to the PPI network
      nodes <- c()
      for (k in seq_len(length(prots))) {
        nodes <- c(nodes, match(prots[k], vertices))
      }
      # Get a subnetwork using the matched proteins
      nodes <- nodes[!is.na(nodes)]
      subnet <- induced_subgraph(net_all, nodes)
      # If TRUE, exclude isolated nodes.
      if (flag_simplify) {
        isolated_nodes <- V(subnet)[degree(subnet)==0]
        subnet_simplified <- delete.vertices(subnet, isolated_nodes)
        subnet <- subnet_simplified
      }
      if (length(V(subnet)) > 0) {
        # Get condition information of subnetwork nodes
        nodes_subnet <- as.vector(V(subnet)$name)
        node_condition <- c()
        for (k in seq_len(length(nodes_subnet))) {
          p <- as.vector(nodes_subnet[k])
          node_condition <- c(node_condition,
                              as.character(conditions[match(p, prots)]))
        }
        # node list
        list_nodes[[i]] <- nodes_subnet
        # condition list
        list_conditions[[i]] <- node_condition
        # plot title
        plot_title <- paste("#",i," (n=", length(V(subnet))," of ",
                            length(conditions),")", sep='')
        titles <- c(titles, plot_title)

        # plot a network
        edge_color <- color_for_graph[i]
        V(subnet)$condition <- node_condition
        vertex_color <- c()
        for (j in seq_len(length(node_condition))) {
          vertex_color <- c(vertex_color,
                            color_for_nodes[match(node_condition[j],
                                                  uniq_condition)])
        }
        E(subnet)$edge.color <- edge_color
        V(subnet)$vertex.color <- vertex_color
        list_subnets[[i]] <- subnet
        plot(subnet, vertex.label.color="black", vertex.label=NA,
             vertex.label.dist=.6, vertex.label.cex=0.8,
             edge.arrow.size=.4, vertex.size=10, vertex.color=vertex_color,
             edge.color=edge_color, vertex.shape=node_shape, main=plot_title,
             layout=layout_nicely(subnet))
      } else {
        plot_NA()
      }
    } else {
      plot_NA()
    }
  }
  dev.off()
  return(list(All_network=net_all,
              Sub_network=list_subnets,
              Data=superset_string,
              num_clusters_in_row=num_clusters_in_row,
              Nodes=list_nodes,
              Conditions=list_conditions,
              Titles=titles,
              out_dir=na_dir,
              is_stringdb=is_stringdb))
}
