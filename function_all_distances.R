## Author: Elsie Horne
## Date created: 21st Jan 2019

## Description: function for direct calculation of Ward's with Gower measure

## Packages:
library(tidyverse)
library(cluster) # for agnes and daisy
library(arrangements) # for combinations

################################################################################
# my_wards function
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
    names[i] <- str_c(c_rem, collapse = " and ")
  }
  names(merges) <- names
  return(merges)
}
################################################################################

################################## check times for euclidean distance
# example data
set.seed(898)
rows <- sample(1:nrow(iris), size = 10)
x <- iris[rows, 1:4] # sample 5 lines of iris as example data
x <- as.data.frame(scale(x)) # standardise to mean = 0 sd = 1
rownames(x) <- as.character(1:nrow(x))

system.time(merges_euc_10 <- my_wards(x, dist = "euclidean"))
# 0.133 seconds
save(merges_euc_10, file = "merges_euc_10.RData")

set.seed(898)
rows <- sample(1:nrow(iris), size = 100)
x <- iris[rows, 1:4] # sample 5 lines of iris as example data
x <- as.data.frame(scale(x)) # standardise to mean = 0 sd = 1
rownames(x) <- as.character(1:nrow(x))

system.time(merges_euc_100 <- my_wards(x, dist = "euclidean"))
# 2403.233 seconds
save(merges_euc_100, file = "merges_euc_100.RData")

################################## check times

# compare gower with my_wards and agnes
set.seed(898)
rows <- sample(1:nrow(iris), size = 10)
x <- iris[rows, 1:4] # sample 5 lines of iris as example data
x <- as.data.frame(scale(x)) # standardise to mean = 0 sd = 1
rownames(x) <- as.character(1:nrow(x))

time_10 <- system.time(merges_gow_10 <- my_wards(x, dist = "gower"))
save(merges_gow_10, file = "merges_gow_10.RData")

dist_gow_10 <- daisy(x, metric = "gower", stand = FALSE)
ag_gow_10 <- agnes(dist_gow_10, method = "ward", stand = FALSE)
sort(ag_gow_10$height)
#plot(ag_gow)

# compare gower with my_wards and agnes
set.seed(898)
rows <- sample(1:nrow(iris), size = 100)
x <- iris[rows, 1:4] # sample 5 lines of iris as example data
x <- as.data.frame(scale(x)) # standardise to mean = 0 sd = 1
rownames(x) <- as.character(1:nrow(x))

time_100 <- system.time(merges_gow_100 <- my_wards(x, dist = "gower"))
save(merges_gow_100, file = "merges_gow_100.RData")

dist_gow_100 <- daisy(x, metric = "gower", stand = FALSE)
ag_gow_100 <- agnes(dist_gow_100, method = "ward", stand = FALSE)
sort(ag_gow_100$height)
#plot(ag_gow)

##################################
clus2 <- names(merges_gow_100[99])
clus2 <- str_split(clus2, " and ")
clus2 <- clus2[[1]]
clus2[1] <- str_split(clus2[1], ",")
clus2[2] <- str_split(clus2[2], ",")
clus2_data1 <- data.frame(ID = clus2[1], clus = 1)
clus2_data2 <- data.frame(ID = clus2[2], clus = 2)
names(clus2_data1) <- c("ID", "clus")
names(clus2_data2) <- c("ID", "clus")
clus2_data <- rbind(clus2_data1, clus2_data2)
clus2_data %>% group_by(clus) %>% count()

##################################

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

