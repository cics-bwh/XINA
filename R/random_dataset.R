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
