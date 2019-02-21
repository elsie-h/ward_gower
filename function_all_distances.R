## Author: Elsie Horne
## Date created: 21st Jan 2019

## Description: function for direct calculation of Ward's with Gower measure

## Packages:
library(tidyverse)
library(cluster) # for agnes and daisy
library(arrangements) # for combinations

# re-write this function to reduce the number of distance matrices that are calculated
# i.e. store all the means of possible merges, and then calculate one distance matrix,
# rather than calculating a new distance matrix for each potential merge

################################################################################
# my_wards function
################################################################################
my_wards <- function(x, dist) {
  # function to compute the error sum of squares for cluster C
  # sum of gower distances between all samples in a cluster to the cluster centroid (mean)
  ess_direct <- function(C) {
    C <- unlist(str_split(C, ","))
    if (length(C) == 1)
      return(0)
    else {
      mean_i <- nrow(x) + 1 # the index for the mean row
      x_mean <- colMeans(x[C, ]) # compute mean of cluster C
      x_C <- rbind(x, x_mean) # samples and mean in one dataset
      d_C <-
        daisy(x_C, metric = "euclidean", stand = FALSE) # compute distances
      d_C <-
        as.matrix(d_C)[mean_i, C] # keep only the row of distances to mean and columns in cluster
      return(sum(d_C * d_C)) # return sum over square of all distances to mean
    }
  }
  # function to compute the error sum of squares for merging two clusters in list L
  change_ess_direct <- function(L) {
    ess_direct(c(L[1], L[2])) - ess_direct(L[1]) - ess_direct(L[2])
  }
  
  levs <- nrow(x) - 1
  merges <- vector(mode = "list", length = levs)
  names <- vector(mode = "list", length = levs)
  clusters <- as.character(1:nrow(x))
  
  for (i in 1:levs) {
    combos <- as.data.frame(t(combinations(x = clusters, k = 2))) %>%
      mutate_all(as.character) # store as character for when clusters get bigger
    d_combos <- lapply(combos, change_ess_direct)
    names(d_combos) <-
      unname(apply(as.matrix(combos), 2, function(x)
        str_c(x, collapse = " and ")))
    d_combos <- unlist(d_combos)
    d_combos <- (2 * d_combos) ^ 0.5 # to match the distances in AGNES
    d_min <- min(d_combos)
    c_rem <- combos[d_combos  == d_min] # clusters to combine
    # merges[i] <- list(d_combos) # store the distance between the merging clusters
    merges[i] <- d_min # if only store the minimum ditance
    c_rem <- as.character(unlist(c_rem))
    c_new <- str_c(unlist(c_rem), collapse = ",")
    clusters <-
      clusters[!(clusters %in% c_rem)] # remove the merged clusters
    clusters <- c(clusters, c_new) # add new merged cluster
    names[i] <- str_c(c_rem, collapse = " and ")
    # merges[i][[1]] <- sort(merges[i][[1]]) # if storing all distances
  }
  names(merges) <- names
  return(merges)
}
################################################################################

# example data
set.seed(898)
rows <- sample(1:nrow(iris), size = 5)
x <- iris[rows, 1:2] # sample 5 lines of iris as example data
x <- as.data.frame(scale(x)) # standardise to mean = 0 sd = 1
rownames(x) <- as.character(1:nrow(x))

merges_gow <- my_wards(x, dist = "gower")
merges_euc <- my_wards(x, dist = "euclidean")

dist_euc <- daisy(x, metric = "euclidean", stand = FALSE)
ag_euc <- agnes(dist_euc, method = "ward", stand = FALSE)
sort(ag_euc$height)
plot(ag_euc)

dist_gow <- daisy(x, metric = "gower", stand = FALSE)
ag_gow <- agnes(dist_gow, method = "ward", stand = FALSE)
sort(ag_gow$height)
plot(ag_gow)

# visualise example data
x %>% ggplot(aes(x = Sepal.Length, y = Sepal.Width, label = rownames(x))) +
  geom_point(shape=15,color="white",size=6)+geom_text()

