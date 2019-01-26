
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

