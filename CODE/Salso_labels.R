library(salso)

# here i upload the matrix of cluster_allocs found in python
clus_allocs <- as.matrix(cluster_allocs)

salso(clus_allocs,loss=binder())

labels=salso(clus_allocs,loss=binder())

labels

