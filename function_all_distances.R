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
    #mean_row <- str_c("mean(", str_c(C, collapse = ","), ")")
    C <- str_c(C, collapse = ",")
    C_ind <- unlist(str_split(C, ","))
    if (length(C_ind) == 1)
      return(0)
    else {
      # keep the clusters as columns and the mean as the row
      d <- d_current[C, C_ind]
      # mean_i <- nrow(x) + 1 # the index for the mean row
      # x_mean <- colMeans(x[C, ]) # compute mean of cluster C
      # x_C <- rbind(x, x_mean) # samples and mean in one dataset
      # d_C <-
      #   daisy(x_C, metric = dist, stand = FALSE) # compute gower distances
      # d_C <-
      #   as.matrix(d_C)[mean_i, C] # keep only the row of distances to mean and columns in cluster
      return(sum(d*d)) # return sum over square of all distances to mean
    }
  }
  # function to compute the error sum of squares for merging two clusters in list L
  change_ess_direct <- function(L) {
    ess_direct(c(L[1], L[2])) - ess_direct(L[1]) - ess_direct(L[2])
  }
  my_mean <- function(combo) {
    combo <- unlist(str_split(combo, ","))
    x_mean <- colMeans(x_current[combo, ])
  }
  
  levs <- nrow(x) - 1
  merges <- vector(mode = "list", length = levs)
  names <- vector(mode = "list", length = levs)
  clusters <- as.character(1:nrow(x))
  x_rows <- as.character(1:nrow(x))
  x_current <- x
  
  for (i in 1:levs) {
    #### calculating change_ess_direct
    # now create just one distance matrix for each level.
    combos <- as.data.frame(t(combinations(x = clusters, k = 2))) %>%
      mutate_all(as.character) # cluster names stored as characters, e.g. "1" or "2,3"
    names(combos) <-
      unname(apply(as.matrix(combos), 2, function(x)
        str_c(x, collapse = ",")))
    # here is where we compute the distance matrix for lev.
    # for each entry in combo, calculate cluster mean and append to x
    means <- sapply(combos, my_mean)
    means <- as.data.frame(t(means))
    x_current <- rbind(x_current, means)
    d_current <- daisy(x_current, metric = dist, stand = FALSE)
    d_current <- as.matrix(d_current)
    # calculate all change_ess_direct for all combos
    d_combos <- lapply(combos, change_ess_direct)
    #### storing results & prepping for next iteration
    names(d_combos) <-
      unname(apply(as.matrix(combos), 2, function(x)
        str_c(x, collapse = " and ")))
    d_combos <- unlist(d_combos)
    d_combos <- (2*d_combos)^0.5 # to match the distances in AGNES
    d_min <- min(d_combos)
    c_rem <- combos[d_combos  == d_min] # clusters to combine
    # merges[i] <- list(d_combos) # store the distance between the merging clusters
    merges[i] <- d_min # if only store the minimum ditance
    c_rem <- as.character(unlist(c_rem))
    c_new <- str_c(unlist(c_rem), collapse = ",")
    clusters <-
      clusters[!(clusters %in% c_rem)] # remove the merged clusters
    clusters <- c(clusters, c_new) # add new merged cluster
    x_rows <- c(x_rows, c_new) # add merged cluster to x_rows
    x_current <- x_current[x_rows,] # new x_current
    names[i] <- str_c(c_rem, collapse = " and ")
    # merges[i][[1]] <- sort(merges[i][[1]]) # if storing all distances
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

set.seed(898)
rows <- sample(1:nrow(iris), size = 100)
x <- iris[rows, 1:4] # sample 5 lines of iris as example data
x <- as.data.frame(scale(x)) # standardise to mean = 0 sd = 1
rownames(x) <- as.character(1:nrow(x))

system.time(merges_euc_100 <- my_wards(x, dist = "euclidean"))
# 2403.233 seconds

################################## check times

# compare gower with my_wards and agnes
set.seed(898)
rows <- sample(1:nrow(iris), size = 10)
x <- iris[rows, 1:4] # sample 5 lines of iris as example data
x <- as.data.frame(scale(x)) # standardise to mean = 0 sd = 1
rownames(x) <- as.character(1:nrow(x))

time_10 <- system.time(merges_gow_10 <- my_wards(x, dist = "gower"))

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

