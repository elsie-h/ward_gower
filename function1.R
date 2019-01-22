## Author: Elsie Horne
## Date created: 21st Jan 2019

## Description: function for direct calculation of Ward's with Gower measure

## Packages:
library(tidyverse)
library(cluster) # for agnes and daisy
library(arrangements) # for combinations

# example data
set.seed(898)
rows <- sample(1:nrow(iris), size = 5)
x <- iris[rows, 1:2] # sample 5 lines of iris as example data
x <- as.data.frame(scale(x)) # standardise to mean = 0 sd = 1
rownames(x) <- as.character(1:nrow(x))
clusters <- as.character(1:nrow(x))

# visualise example data
x %>% ggplot(aes(x = Sepal.Length, y = Sepal.Width, label = rownames(x))) +
  geom_point(shape=15,color="white",size=6)+geom_text()


# function to compute the error sum of squares for cluster C
# sum of gower distances between all samples in a cluster to the cluster centroid (mean)
ess_direct <- function(C) {
  C <- unlist(str_split(C, ","))
  if (length(C) == 1) return(0) else {
    mean_i <- nrow(x) + 1 # the index for the mean row
    x_mean <- colMeans(x[C,]) # compute mean of cluster C
    x_C <- rbind(x, x_mean) # samples and mean in one dataset
    d_C <- daisy(x_C, metric = "euclidean", stand = FALSE) # compute gower distances
    d_C <- as.matrix(d_C)[mean_i,C] # keep only the row of distances to mean and columns in cluster
    return(sum(d_C*d_C)) # return sum over square of all distances to mean
  }
}
# function to compute the error sum of squares for merging two clusters in list L
change_ess_direct <- function(L) {
  ess_direct(c(L[1], L[2])) - ess_direct(L[1]) - ess_direct(L[2])
}

levs <- nrow(x) - 1
merges <- vector(mode = "list", length = levs)
names <- vector(mode = "list", length = levs)

for (i in 1:levs) {
  combos <- as.data.frame(t(combinations(x = clusters, k = 2))) %>%
    mutate_all(as.character) # store as character for when clusters get bigger
  d_combos <- lapply(combos, change_ess_direct)
  d_min <- min(unlist(d_combos))
  c_rem <- combos[d_combos  == d_min] # clusters to combine
  merges[i] <-  d_min # store the distance between the merging clusters
  c_rem <- as.character(unlist(c_rem))
  c_new <- str_c(unlist(c_rem), collapse = ",")
  clusters <- clusters[!(clusters %in% c_rem)] # remove the merged clusters
  clusters <- c(clusters, c_new) # add new merged cluster
  names[i] <- str_c(c_rem, collapse = " and ")
}
names(merges) <- names