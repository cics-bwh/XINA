% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering_plots.R
\name{mutate_colors}
\alias{mutate_colors}
\title{mutate_colors}
\usage{
mutate_colors(condition_composition, color_for_condition,
  null_color = "gray", threshold_percent = 50)
}
\arguments{
\item{condition_composition}{A data frame generated by \link[XINA]{plot_condition_compositions}}

\item{color_for_condition}{A vector like 'color_for_condition' of \link[XINA]{xina_clustering}}

\item{null_color}{Default is 'gray'. This color is for clusters that are not biased to any of experimental conditions}

\item{threshold_percent}{Default is 50.  The percentage threshold for giving new colors}
}
\value{
A data frame containing statistics of XINA network centrality scores
}
\description{
'mutate_colors' generates new color scheme for XINA clustering plot based on condition composition results (\link[XINA]{plot_condition_compositions}).
If any clusters have higher percentage than the 'threshold_percent', XINA will assign new colors in accordance to 'color_for_condition'.
If not, XINA will give 'gray' color or user-defined color via 'null_color' parameter.
}
\examples{
# load XINA example data
data(xina_example)

# Plot condition composition pie-chart with default option
condition_composition <- plot_condition_compositions(example_clusters)
example_clusters$color_for_clusters <- mutate_colors(condition_composition,
example_clusters$color_for_condition)
plot_clusters(example_clusters, xval=c(0,2,6,12,24,48,72), xylab=FALSE)

}
