<<<<<<< HEAD
##############################################
# The deveroper and the maintainer:          #
#   Lang Ho Lee (lhlee@bwh.harvard.edu)      #
#   Sasha A. Singh (sasingh@bwh.harvard.edu) #
##############################################
=======
#############################################
# The deveroper and the maintainer:         #
#   Lang Ho Lee (lhlee@bwh.harvard.edu)     #
#   Sasha A. Singh (sasingh@bwh.harvard.edu)#
#############################################
>>>>>>> 897b60ba80ec40f31146fc7e7fffb563da4f3b87

# Mute warnings
options(warn=1)

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
<<<<<<< HEAD
#' @param fill_color Default is 'darkgray'. You can change color of bars.
=======
>>>>>>> 897b60ba80ec40f31146fc7e7fffb563da4f3b87
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
<<<<<<< HEAD
                                    sig_score="pvalue", num_terms=0, get_log=TRUE,
                                    fill_color='darkgray'){
=======
                                    sig_score="pvalue", num_terms=0, get_log=TRUE){
>>>>>>> 897b60ba80ec40f31146fc7e7fffb563da4f3b87
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
<<<<<<< HEAD
    ggplot(aes(term_description, Enriched_score)) + geom_bar(stat="identity", fill=fill_color) + coord_flip() + ylab(y_lab)
=======
    ggplot(aes(term_description, Enriched_score)) + geom_bar(stat="identity") + coord_flip() + ylab(y_lab)
>>>>>>> 897b60ba80ec40f31146fc7e7fffb563da4f3b87
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
