# ward_gower

I've come across several papers which use Gower's distance with Ward's hierarchical clustering due to mixed data types, implemented by cluster::agnes in R. Ward's algorithm is intended for use with the squared Euclidean distance, and cluster::agnes makes use of some properties of the Euclidean distance for computational efficiency. 'Writeup.Rmd' demonstrates that Ward's method with the Gower distance is not implemented properly with cluster::agnes in R.
