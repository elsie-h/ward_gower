## Author: Elsie Horne
## Date created: 21st Jan 2019

## Description: function for direct calculation of Ward's

## Packages:
require(tidyverse)
require(cluster) # for daisy
require(arrangements) # for combinations

################################################################################
# my_wards function:
# wards clustering by calculating centres and distances directly
################################################################################
my_wards <- function(x, dist) {
  # function to compute the ESS for cluster C
  # sum of distances between all samples in a cluster to the cluster centroid
  ess_direct <- function(C) {
    C <- str_c(C, collapse = ",") # samples in cluster together e.g. "1,2,3"
    C_ind <- unlist(str_split(C, ",")) # samples in cluster individually e.g. "1" "2" "3"
    if (length(C_ind) == 1) # if singleton cluster ESS = 0
      return(0)
    else {
      # keep the samples as columns and cluster centroid as row
      # calculate distance matrix
      d <- d_current[C, C_ind]
      return(sum(d*d)) # sum over square of all distances to mean
    }
  }
  # function to compute the change in ESS for merging two clusters in list L
  change_ess_direct <- function(L) {
    ess_direct(c(L[1], L[2])) - ess_direct(L[1]) - ess_direct(L[2])
  }
  # function to take column means for all possible cluster merges
  my_mean <- function(combo) {
    combo <- unlist(str_split(combo, ","))
    x_mean <- colMeans(x_current[combo, ])
  }
  
  levs <- nrow(x) - 1 # number of levels in hierarchy
  merges <- vector(mode = "list", length = levs) # output of below loop
  names <- vector(mode = "list", length = levs) # names for output
  clusters <- as.character(1:nrow(x)) # current clusters (singletons at start)
  x_rows <- as.character(1:nrow(x)) # rows of x
  x_current <- x # data matrix for current level
  
  for (i in 1:levs) {
    # all possible merges
    combos <- as.data.frame(t(combinations(x = clusters, k = 2))) %>%
      mutate_all(as.character) # cluster names stored as characters, e.g. "2,3"
    names(combos) <-
      unname(apply(as.matrix(combos), 2, function(x)
        str_c(x, collapse = ",")))
    # for each entry in combo, calculate cluster mean and append to x
    means <- sapply(combos, my_mean)
    means <- as.data.frame(t(means))
    x_current <- rbind(x_current, means)
    # compute distance matrix for this level
    d_current <- daisy(x_current, metric = dist, stand = FALSE)
    d_current <- as.matrix(d_current)
    # calculate all change_ess_direct for all combos
    d_combos <- lapply(combos, change_ess_direct)
    # storing results & prep for next iteration
    names(d_combos) <-
      unname(apply(as.matrix(combos), 2, function(x)
        str_c(x, collapse = " and ")))
    d_combos <- unlist(d_combos)
    d_combos <- (2*d_combos)^0.5 # to match the distances in AGNES
    d_min <- min(d_combos)
    c_rem <- combos[d_combos  == d_min] # clusters to combine
    merges[i] <- d_min # if only store the minimum ditance
    c_rem <- as.character(unlist(c_rem))
    c_new <- str_c(unlist(c_rem), collapse = ",")
    clusters <-
      clusters[!(clusters %in% c_rem)] # remove the merged clusters
    clusters <- c(clusters, c_new) # add new merged cluster
    x_rows <- c(x_rows, c_new) # add merged cluster to x_rows
    x_current <- x_current[x_rows,] # new x_current
    #names[i] <- str_c(c_rem, collapse = " and ")
    names[i] <- str_c(clusters, collapse = "--")
  }
  names(merges) <- names
  return(merges)
}
################################################################################


################################################################################
# results
# use output from m_wards to calculate the cluster solution for a given k
################################################################################
my_clusters <- function(my_wards_out, k) {
  r <- names(my_wards_out)
  r <- r[[nrow(x) - k]][[1]]
  r <- str_split(r, "--")[[1]]
  r <- sapply(r, function(x) str_split(x, ","))
  r <- sapply(r, as.numeric)
  
  x_r <- data.frame()
  for (i in 1:length(r)) {
    m_r <- matrix(data = c(r[[i]], rep(i, length = length(r[[i]]))), ncol = 2)
    m_r <- as.data.frame(m_r)
    x_r <- rbind(x_r, m_r)
  }
  x_r <- unname(unlist(arrange(x_r, V1)[2]))
  return(x_r)
}
################################################################################
