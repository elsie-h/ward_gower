# Introduction to statistical learning example
# 10.6 Lab 3: NCI60 Data Example

library(ISLR)
library(cluster)

nci.labs=NCI60$labs
nci.data=NCI60$data
dim(nci.data)

### PCA on NCI60
# first scale
pr.out=prcomp(nci.data, scale=TRUE)


# function to assign colour to each cell line based on corresponding cancer type
Cols=function(vec){
  cols=rainbow(length(unique(vec)))
  return(cols[as.numeric(as.factor(vec))])
}

# Plot first 3 component score vectors
par(mfrow=c(1,2))
plot(pr.out$x[,1:2], col=Cols(nci.labs), pch=19,
       xlab="Z1",ylab="Z2")
plot(pr.out$x[,c(1,3)], col=Cols(nci.labs), pch=19,
       xlab="Z1",ylab="Z3")

summary(pr.out)

# scree plot
pve=100*pr.out$sdev^2/sum(pr.out$sdev^2)
par(mfrow=c(1,2))
plot(pve, type="o", ylab="PVE", xlab="Principal Component",
       col =" blue ")
plot(cumsum(pve), type="o", ylab="Cumulative PVE", xlab="
       Principal Component ", col =" brown3 ")

# the first 7 components explain about 40% of total of vaiance
# scree plot suggests after 7 there is a marked decrease in the 
# proportion of variance explained.

###
# hierarchical clustering on first 5 component scores
hc.out=hclust(dist(pr.out$x[,1:5]))
plot(hc.out, labels=nci.labs, main="Hier. Clust. on First
       Five Score Vectors ")
table(cutree(hc.out,4), nci.labs)

# with agnes


x<- as.data.frame(pr.out$x[,1:5])
rownames(x) <- as.character(1:nrow(x))
x_ <- x %>%
  mutate(x_row = as.numeric(row.names(x))) %>%
  full_join(x_r, by = "x_row")


ag_euc <- agnes(x, method = "ward", stand = FALSE)
sort(ag_euc$height)
plot(ag_euc)
table(cutree(ag_euc,4), nci.labs)

time_euc <- system.time(mw_euc <- my_wards(x = x, dist = "euclidean"))
mw_euc

# results are same with my_wards and agnes


d <- daisy(x, metric = "gower", stand = FALSE)
ag_gow <- agnes(d, method = "ward", stand = FALSE)
sort(ag_gow$height)
plot(ag_gow)
table(cutree(ag_gow,5), nci.labs)

time_gow <- system.time(mw_gow <- my_wards(x = x, dist = "gower"))
unname(unlist(mw_gow))
my_5 <- my_clusters(mw_gow, k = 5)
agnes_5 <- cutree(ag_gow,5)
table(my_5, agnes_5)

plot_data <- data.frame(agnes = sort(ag_gow$height), my = unname(unlist(mw_gow)))

sort(ag_gow$height)/unname(unlist(mw_gow))

k <- 5
names(mw_gow[nrow(x) + 1 - k])
