## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----installation, eval=FALSE--------------------------------------------
#  # XINA requires R >= 3.5.1
#  install.packages('devtools')
#  library('devtools')
#  install_github('langholee/XINA')

## ----import libraries----------------------------------------------------
library(XINA)
library(igraph)
library(ggplot2)
library(STRINGdb)

## ----example random dataset----------------------------------------------
# Generate random multiplexed time-series data
random_data_info <- make_random_xina_data()

# The number of proteins
random_data_info$size

# Time points
random_data_info$time_points

# Three conditions
random_data_info$conditions

## ----check randomly generated data files---------------------------------
Control <- read.csv("Control.csv", check.names=FALSE, stringsAsFactors = FALSE)
Stimulus1 <- read.csv("Stimulus1.csv", check.names=FALSE, stringsAsFactors = FALSE)
Stimulus2 <- read.csv("Stimulus2.csv", check.names=FALSE, stringsAsFactors = FALSE)

head(Control)
head(Stimulus1)
head(Stimulus2)

## ----data matrix---------------------------------------------------------
head(Control[random_data_info$time_points])

## ----fix random seed-----------------------------------------------------
set.seed(0)

## ----set up for the clustering-------------------------------------------
# Data files
data_files <- paste(random_data_info$conditions, ".csv", sep='')
data_files

# time points of the data matrix
data_column <- random_data_info$time_points
data_column

## ----XINA clustering-----------------------------------------------------
# Run the model-based clusteirng
clustering_result <- xina_clustering(data_files, data_column=data_column, nClusters=30)

## ----XINA clustering with VVI covariance matrix, eval=FALSE--------------
#  # Model-based clustering using VVI covariance model
#  clustering_result_vvi <- xina_clustering(data_files, data_column=data_column, nClusters=30, chosen_model='VVI')

## ----XINA with k-means clustering----------------------------------------
clustering_result_km <- xina_clustering(data_files, data_column=data_column, nClusters=30, chosen_model='kmeans')

## ----load previous data--------------------------------------------------
clustering_result_reloaded <- load_previous_results(".")
head(clustering_result_reloaded$aligned)

## ----load xina_example---------------------------------------------------
data(xina_example)

## ----XINA clustering plot------------------------------------------------
plot_clusters(example_clusters, xylab=FALSE)

## ----clustering plot2, eval=FALSE----------------------------------------
#  plot_clusters(example_clusters, xval=c(0,2,6,12,24,48,72), xylab=FALSE)
#  
#  # You can change Y axis limits to be ranged from 0 to 0.35
#  plot_clusters(example_clusters, y_lim=c(0,0.35), xval=c(0,2,6,12,24,48,72), xylab=FALSE)
#  
#  # If you need, you can modify the clustering plot by creating your own ggplot theme.
#  theme1 <- theme(title=element_text(size=8, face='bold'),
#                  axis.text.x = element_text(size=7),
#                  axis.text.y = element_blank(),
#                  axis.ticks.x = element_blank(),
#                  axis.ticks.y = element_blank(),
#                  axis.title.x = element_blank(),
#                  axis.title.y = element_blank())
#  plot_clusters(example_clusters, ggplot_theme=theme1)

## ----XINA condition composition------------------------------------------
condition_composition <- plot_condition_compositions(example_clusters)
tail(condition_composition)

## ----XINA condition composition2, eval=FALSE-----------------------------
#  theme2 <- theme(legend.position="none", title=element_text(size=7, face='bold'))
#  condition_composition <- plot_condition_compositions(example_clusters, ggplot_theme=theme2)
#  
#  # Or you can adjust the legend size.
#  theme3 <- theme(legend.key.size = unit(0.3, "cm"),
#                  legend.text=element_text(size=5),
#                  title=element_text(size=7, face='bold'))
#  condition_composition <- plot_condition_compositions(example_clusters, ggplot_theme=theme3)
#  
#  # You can utilize your own colors for condition composition charts.
#  # Make a new color code for conditions
#  condition_colors <- c("tomato","steelblue1","gold")
#  names(condition_colors) <- example_clusters$condition
#  example_clusters$color_for_condition <- condition_colors
#  
#  # Draw condition composition pie-chart with the new color code
#  condition_composition <- plot_condition_compositions(example_clusters, ggplot_theme=theme3)
#  
#  # XINA also can draw bull's-eye plots instead of pie-charts.
#  condition_composition <- plot_condition_compositions(example_clusters, ggplot_theme=theme3, bullseye = TRUE)

## ----clustering plot with condition colors-------------------------------
# New colors for the clustering plot based on condition composition
example_clusters$color_for_clusters <- mutate_colors(condition_composition, example_clusters$color_for_condition)
example_clusters$color_for_clusters

## ----drawing clustering plot with the mutated colors, eval=FALSE---------
#  plot_clusters(example_clusters, xval=c(0,2,6,12,24,48,72), xylab=FALSE)

## ----clustering plot with lowered threshold, eval=FALSE------------------
#  example_clusters$color_for_clusters <- mutate_colors(condition_composition, example_clusters$color_for_condition, threshold_percent=40)
#  plot_clusters(example_clusters, xval=c(0,2,6,12,24,48,72), xylab=FALSE)

## ----XINA comigration search---------------------------------------------
classes <- as.vector(example_clusters$condition)
classes

all_cor <- alluvial_enriched(example_clusters, classes)
head(all_cor)

## ----XINA comigration search with different order, eval=FALSE------------
#  all_cor_Stimulus1_start <- alluvial_enriched(example_clusters, c(classes[2],classes[1],classes[3]))
#  head(all_cor_Stimulus1_start)

## ----XINA comigration search with comigration size filter, eval=FALSE----
#  cor_bigger_than_10 <- alluvial_enriched(example_clusters, classes, comigration_size=10)
#  head(cor_bigger_than_10)

## ----finding the condition-biased clusters, eval=FALSE-------------------
#  condition_biased_comigrations <- get_condition_biased_comigrations(clustering_result=example_clusters,
#                                                                     count_table=cor_bigger_than_10,
#                                                                     selected_conditions=classes,
#                                                                     condition_composition=condition_composition,
#                                                                     threshold_percent=50, color_for_null='gray',
#                                                                     color_for_highly_matched='yellow', cex = 1)

## ----protein list of specific comigrations-------------------------------
condition_biased_clusters <- condition_composition[condition_composition$Percent_ratio>=50,]
control_biased_clusters <- condition_biased_clusters[condition_biased_clusters$Condition=="Control",]$Cluster
Stimulus1_biased_clusters <- condition_biased_clusters[condition_biased_clusters$Condition=="Stimulus1",]$Cluster
Stimulus2_biased_clusters <- condition_biased_clusters[condition_biased_clusters$Condition=="Stimulus2",]$Cluster

# Get the proteins showing condition-specific expression patterns in three conditions
proteins_found_in_condition_biased_clusters <- subset(example_clusters$aligned, Control==control_biased_clusters[1] & Stimulus1==Stimulus1_biased_clusters[1] & Stimulus2==Stimulus2_biased_clusters[1])
nrow(proteins_found_in_condition_biased_clusters)
head(proteins_found_in_condition_biased_clusters)

## ----comigrations with comigration size and significance filters---------
Control_Stimulus1_significant <- alluvial_enriched(example_clusters, c("Control","Stimulus1"), comigration_size=5, pval_threshold=0.05)
head(Control_Stimulus1_significant)

## ----STRING DB set up----------------------------------------------------
string_db <- STRINGdb$new( version="10", species=9606, score_threshold=0, input_directory="" )
string_db

## ----select PPI confidence level of STRING DB, eval=FALSE----------------
#  get.edge.attribute(sub_graph)$combined_score
#  
#  # Medium confidence PPIs only
#  string_db_med_confidence <- subgraph.edges(string_db$graph, which(E(string_db$graph)$combined_score>=400), delete.vertices = TRUE)
#  
#  # High confidence PPIs only
#  string_db_high_confidence <- subgraph.edges(string_db$graph, which(E(string_db$graph)$combined_score>=700), delete.vertices = TRUE)

## ----XINA analysis with STRING DB----------------------------------------
xina_result <- xina_analysis(example_clusters, string_db)

## ----HPRD PPI network----------------------------------------------------
# Construct HPRD PPI network
data(HPRD)
ppi_db_hprd <- simplify(graph_from_data_frame(d=hprd_ppi, directed=FALSE), remove.multiple = FALSE, remove.loops = TRUE)
head(hprd_ppi)

## ----XINA analysis with HPRD---------------------------------------------
# Run XINA with HPRD protein-protein interaction database
xina_result_hprd <- xina_analysis(example_clusters, ppi_db_hprd, is_stringdb=FALSE)

## ----plotting PPI networks of all the clusters, eval=FALSE---------------
#  # XINA network plots labeled gene names
#  xina_plots(xina_result, example_clusters)

## ----xina_plots without labels, eval=FALSE-------------------------------
#  # XINA network plots without labels
#  xina_plots(xina_result, example_clusters, vertex_label_flag=FALSE)

## ----xina_plots with tree layout, eval=FALSE-----------------------------
#  # Plot PPI networks with tree layout
#  xina_plots(xina_result, example_clusters, vertex_label_flag=FALSE, layout_specified='tree')

## ----xina_plot_bycluster-------------------------------------------------
xina_plot_bycluster(xina_result, example_clusters, cl=1)

## ----xina_plot_bycluster2, eval=FALSE------------------------------------
#  xina_plot_bycluster(xina_result, example_clusters, cl=1, condition="Control")

## ----xina_plot_all-------------------------------------------------------
img_size <- c(3000,3000)
xina_plot_all(xina_result, example_clusters, img_size=img_size)

## ----xina_plot_all2, eval=FALSE------------------------------------------
#  xina_plot_all(xina_result, example_clusters, condition="Control", img_size=img_size)
#  xina_plot_all(xina_result, example_clusters, condition="Stimulus1", img_size=img_size)
#  xina_plot_all(xina_result, example_clusters, condition="Stimulus2", img_size=img_size)

## ----xina_plot_all_eigen, eval=FALSE-------------------------------------
#  xina_plot_all(xina_result, example_clusters, centrality_type="Eigenvector",
#                edge.color = 'gray', img_size=img_size)
#  xina_plot_all(xina_result, example_clusters, condition="Control", centrality_type="Eigenvector",
#                edge.color = 'gray', img_size=img_size)
#  xina_plot_all(xina_result, example_clusters, condition="Stimulus1", centrality_type="Eigenvector",
#                edge.color = 'gray', img_size=img_size)
#  xina_plot_all(xina_result, example_clusters, condition="Stimulus2", centrality_type="Eigenvector",
#                edge.color = 'gray', img_size=img_size)

## ----enrichment tests, eval=FALSE----------------------------------------
#  enriched_fdr_0.05 <- xina_enrichment(xina_result, example_clusters, string_db, pval_threshold=0.05)

## ----network analysis of comigrated proteins-----------------------------
protein_list <- proteins_found_in_condition_biased_clusters$`Gene name`
protein_list

plot_title_ppi <- "PPI network of the proteins found in the condition-biased clusters"

# Draw PPI networks and compute eigenvector centrality score.
comigrations_condition_biased_clusters <- xina_plot_single(xina_result, protein_list, centrality_type="Eigenvector",
                                                           vertex.label.cex=0.5, vertex.size=8, vertex.label.dist=1,
                                                           main=plot_title_ppi)

## ----three breaks, results="hide", eval=FALSE----------------------------
#  # Draw with less breaks
#  comigrations_condition_biased_clusters <- xina_plot_single(xina_result, protein_list, centrality_type="Eigenvector",
#                                                             vertex.label.cex=0.5, vertex.size=8, vertex.label.dist=1,
#                                                             num_breaks=3, main=plot_title_ppi)
#  

## ----test of xina network plotting, eval=FALSE---------------------------
#  # Draw without labels
#  comigrations_condition_biased_clusters <- xina_plot_single(xina_result, protein_list, centrality_type="Eigenvector",
#                                                             vertex.size=10, vertex_label_flag=FALSE)

## ----layout test, eval=FALSE---------------------------------------------
#  # Sphere layout
#  comigrations_condition_biased_clusters <- xina_plot_single(xina_result, protein_list, centrality_type="Eigenvector",
#                                                             vertex.label.cex=0.5, vertex.size=8, vertex.label.dist=1,
#                                                             main=plot_title_ppi, layout_specified = "sphere")
#  # Star layout
#  comigrations_condition_biased_clusters <- xina_plot_single(xina_result, protein_list, centrality_type="Eigenvector",
#                                                             vertex.label.cex=0.5, vertex.size=8, vertex.label.dist=1,
#                                                             main=plot_title_ppi, layout_specified = "star")
#  
#  # Gem layout
#  comigrations_condition_biased_clusters <- xina_plot_single(xina_result, protein_list, centrality_type="Eigenvector",
#                                                             vertex.label.cex=0.5, vertex.size=8, vertex.label.dist=1,
#                                                             main=plot_title_ppi, layout_specified = "gem")
#  
#  # Tree layout
#  comigrations_condition_biased_clusters <- xina_plot_single(xina_result, protein_list, centrality_type="Eigenvector",
#                                                             vertex.label.cex=0.5, vertex.size=8, vertex.label.dist=1,
#                                                             main=plot_title_ppi, layout_specified = "tree")
#  
#  # Circle layout
#  comigrations_condition_biased_clusters <- xina_plot_single(xina_result, protein_list, centrality_type="Eigenvector",
#                                                             vertex.label.cex=0.5, vertex.size=8, vertex.label.dist=1,
#                                                             main=plot_title_ppi, layout_specified = "circle")
#  
#  # Random layout
#  comigrations_condition_biased_clusters <- xina_plot_single(xina_result, protein_list, centrality_type="Eigenvector",
#                                                             vertex.label.cex=0.5, vertex.size=8, vertex.label.dist=1,
#                                                             main=plot_title_ppi, layout_specified = "random")
#  
#  # Nicely layout
#  comigrations_condition_biased_clusters <- xina_plot_single(xina_result, protein_list, centrality_type="Eigenvector",
#                                                             vertex.label.cex=0.5, vertex.size=8, vertex.label.dist=1,
#                                                             main=plot_title_ppi, layout_specified = "nicely")
#  

## ----enrichment test of the comigrated proteins--------------------------
# Enrichment test using KEGG pathway terms
kegg_enriched <- xina_enrichment(string_db, protein_list, enrichment_type = "KEGG", pval_threshold=0.1)
head(kegg_enriched$KEGG)
# plot enrichment test results
plot_enrichment_results(kegg_enriched$KEGG, num_terms=20)

## ----enrichment test of the comigrated proteins2, eval=FALSE-------------
#  # Enrichment test using GO pathway terms
#  go_enriched <- xina_enrichment(string_db, protein_list, enrichment_type = "GO", pval_threshold=0.1)
#  head(go_enriched$Process)
#  head(go_enriched$Function)
#  head(go_enriched$Component)

## ----draw GO enrichment results, eval=FALSE------------------------------
#  plot_enrichment_results(go_enriched$Process, num_terms=20)
#  plot_enrichment_results(go_enriched$Function, num_terms=20)
#  plot_enrichment_results(go_enriched$Component, num_terms=20)

