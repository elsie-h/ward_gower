## Author: Elsie Horne
## Date created: 16th Jan 2019

## Description: rough work for the 'Gower_proof.Rmd'

## Packages:
library(tidyverse)
library(cluster)
library(arrangements)


set.seed(898)
rows <- sample(1:nrow(iris), size = 5)
x <- iris[rows, 1:4] # sample 5 lines of iris as example data
d <- daisy(x, metric = "gower", stand = TRUE)

# First merge is based on smallest distance between samples.
# Merge samples 90 (row 4) & 70 (row 3)
# Create new line for x which is arithmetic mean of samples 90 & 70
# Let cluster A = sample 90, cluster B = sample 70, 
# cluster R = merged clusters A and B, cluster Q = sample 43
# |A| = |B| = |Q| = 1, |R| = 2

xm <- colMeans(x[c(3,4),]) # x merged

xn <- rbind(x[-c(3,4),], xm) # x new

dn <- daisy(xn, metric = "gower", stand = TRUE)

d2RQ <- (2*2*1/(2+1))*0.6346158^2 # calculated directly from new distance matrix
d2AQ <- (2*1*1/(1+1))*0.64436147^2
d2BQ <- (2*1*1/(1+1))*0.62487013^2
d2AB <- (2*1*1/(1+1))*0.03949134^2

test <- ((1+1)/(2+1))*d2AQ + ((1+1)/(2+1))*d2BQ - (1/(2+1))*d2AB # calculated using update equation

d2RQ == test

# quite close, I should write a function for doing this to check over multiple merges.
##################################################


x %>% ggplot(aes(x = Sepal.Length, y = Sepal.Width, label = rownames(x))) +
  geom_point(shape=15,color="white",size=6)+geom_text()

dg <- daisy(x, metric = "euclidean", stand = TRUE)
ag <- agnes(x, method = "ward", stand = TRUE)
sort(ag$height)
plot(ag)



#ward_gower <- function(x) {
  #### first merge
  




  #for (i in 1:levs) {

    ess_direct <- function(C) {
      C <- unlist(str_split(C, ","))
      if (length(C) == 1) return(0) else {
      mean_i <- nrow(x) + 1 # the index for the mean row
      x_mean <- colMeans(x[C,]) # compute mean of cluster C
      x_C <- rbind(x, x_mean) # samples and mean in one dataset
      d_C <- daisy(x_C, metric = "gower", stand = FALSE) # compute gower distances
      d_C <- as.matrix(d_C)[mean_i,C] # keep only the row of distances to mean and columns in cluster
      return(sum(d_C*d_C)) # return sum over square of all distances to mean
      }
    }
    change_ess_direct <- function(C) {
      ess_direct(c(C[1], C[2])) - ess_direct(C[1]) - ess_direct(C[2])
    }
    
    
    set.seed(898)
    rows <- sample(1:nrow(iris), size = 5)
    x <- iris[rows, 1:2] # sample 5 lines of iris as example data
    x <- as.data.frame(scale(x)) # standardise to mean = 0 sd = 1
    rownames(x) <- as.character(1:nrow(x))
    clusters <- as.character(1:nrow(x))
    
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
    

  #}
  
#############################################  
  d <- daisy(x, metric = "gower", stand = TRUE)
  dmin <- min(d)
  d <- as.matrix(d)
  d <- d*lower.tri(d)
  dd <- as.data.frame(d)
  dd <- cbind(rows, dd)
  dd <- filter_all(dd, any_vars(. == dmin))
  cA <- as.character(unname(unlist(select(dd, rows))))
  cB <- as.character(rows[unlist(dd[-1]) == dmin])
  
  # x merged
  xmerge <- matrix(colMeans(x[c(cA,cB),]), 
                   nrow = 1, 
                   dimnames = list("R", names(x))) 
  # x new
  xnew <- rbind(x[-which(rownames(x) %in% c(cA, cB)),], xm) # x new
  
  
  for (i in 1:nrow(x)) {
    
    # last line
    dold <- dnew
  }
}
test <- vector("list", 5)
for (i in 1:5) {
  print("hello")
  print(data.frame(one  = i, two = i))
  test[i] <-list(one = c(i, i))
}


set.seed(898)
rows <- sample(1:nrow(iris), size = 5)
x <- iris[rows, 1:4] # sample 5 lines of iris as example data


d <- daisy(x, metric = "gower", stand = TRUE)
dmin <- min(d)
d <- as.matrix(d)
d <- d*lower.tri(d)
dd <- as.data.frame(d)
dd <- cbind(rows, dd)
dd <- filter_all(dd, any_vars(. == dmin))
cA <- as.character(unname(unlist(select(dd, rows))))
cB <- as.character(rows[unlist(dd[-1]) == dmin])

xm <- matrix(colMeans(x[c(cA,cB),]), nrow = 1, dimnames = list("R", names(x))) # x merged
xn <- rbind(x[-which(rownames(x) %in% c(cA, cB)),], xm) # x new
dn <- daisy(xn, metric = "gower", stand = TRUE) # new dissimilarity matrix
dn <- as.matrix(dn)
dn <- dn*lower.tri(dn)

nA <- 1
nB <- 1
nR <- 2
nQ <- 1

cQ <- "67" # cluster for next merge

# values required for update equation:

dcc <- function(c1, c2) {
  args <- match.call()
  args <- c(args[2][[1]], args[3][[1]])
  if (any(args %in% c("cA", "cB"))) dist <- d else dist <- dn
  dcc <- dist[c1, c2]
  if (identical(dcc, 0)) dcc <- dist[c2, c1]
  return(dcc)
}

dRQ <- dcc("R", cQ)
dAB <- dcc(cA, cB)
dAQ <- dcc(cA, cQ)
dBQ <- dcc(cB, cQ)

# LHS
d2RQ <- (2*nR*nQ/(nR+nQ))*dRQ^2 # calculated directly from new distance matrix
# RHS
d2AQ <- (2*nA*nQ/(nA+nQ))*dAQ^2
d2BQ <- (2*nB*nQ/(nB+nQ))*dBQ^2
d2AB <- (2*nA*nB/(nA+nB))*dAB^2

update <- ((nA+nQ)/(nR+nQ))*d2AQ + ((nB+nQ)/(nR+nQ))*d2BQ - (nQ/(nR+nQ))*d2AB # calculated using update equation

update
d2RQ

#### carry on from here

unlist(d)
str(d)
as.vector(d)
dmin <- min(d)
d[1,2]
rownames(d)
