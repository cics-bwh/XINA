#############################################
# The deveroper and the maintainer:         #
#   Lang Ho Lee (lhlee@bwh.harvard.edu)     #
#   Sasha A. Singh (sasingh@bwh.harvard.edu)#
#############################################

# Mute warnings
options(warn=1)

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
