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

(d <- daisy(x, metric = "euclidean"))
ag <- agnes(d, method = "ward")
sort(ag$height) # first merge is 3 and 4

# calculate directly using eq.22 & eq.24 from AGNES
# squared euclidean distance between samples 3 & 4 (because singleton)
(d2_3_4  <- 0.1128665*0.1128665 ) # eq.22
(ess_3_4 <- 0.5*d2_3_4) #eq. 24
# ess_3_4 matches the merge calculated by my function

(2*ess_3_4)^0.5

d2_3_1 <- 2.3859301*2.3859301
d2_4_1 <- 2.3237188*2.3237188

(d_1_34 <- (2/3)*d2_3_1 + (2/3)*d2_4_1 - (1/3)*d2_3_4)

(ess_1_34 <- 0.5*7.390641)
# ess_1_34 also matches the one by my function

# but neither of these match the one calculated by AGNES...

d_1_34^0.5

merge_1 <- unlist(merges[[1]])
(2*merge_1)^0.5


(2*mat[,1])

apply(mat, 2, sum)
function(x) (2*x)^0.5)

