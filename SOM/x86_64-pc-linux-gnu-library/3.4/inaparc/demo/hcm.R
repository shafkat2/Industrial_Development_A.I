# Demo for K-means clustering using kmpp and inofrep prototype initialization
# functions in the package 'inaparc'
#

# Install and/or load the required libraries
for(package in c("inaparc", "kpeaks", "cluster", "factoextra", "GGally", "NbClust")){
  if(!require(package, character.only=TRUE, quietly=TRUE)){
    install.packages(package)
  }
  library(package, character.only=TRUE)
}

# Data loading and pre-processing

# Use Iris data set
data(iris)

# Visualize X using GGaly
GGally::ggpairs(iris, aes(colour = Species))

# Pick the features except the 5th column 
X <- iris[,-5]

# Convert type of X to the matrix class
X <- as.matrix(X)
str(X)
head(X)
tail(X)

# Omit the NAs if any
X <- na.omit(X)

# For standardization, uncomment the following two lines
# X <- scale(X, center = TRUE, scale = TRUE) 
# head(X); tail(X)

# Compute the distance matrix 
dist.x <- factoextra::get_dist(X, stand = TRUE, method = "pearson")

par(ask=TRUE)
# Visualize the distance matrix of X
factoextra::fviz_dist(dist.x, 
   gradient = list(low = "#00AFBB", 
   mid = "white", high = "#FC4E07"))

# Quick estimation of the number of clusters by using the findk function of the kpeaks package
nc1 <- kpeaks::findk(X, nbins=30)$mtlk
cat("The optimal number of clusters (k) suggested by kpeaks is: ", nc1, "\n")

# Or finding the optimal number of clusters by using the internal indexes in the NbClust package
nbclust.res <- NbClust::NbClust(X, 
   distance = "euclidean",
   min.nc = 2, max.nc = 10, 
   method = "ward.D", 
   index ="all")

# Finding the most repeated number of clusters and assign it as the optimal number of clusters 
tabnc <- table(nbclust.res$Best.nc[1,])
nc2 <- as.integer(names(tabnc)[which.max(tabnc)])
cat("The optimal number of clusters (k) suggested by NbClust is: ", nc2, "\n")

par(mfrow=c(1,1))
# Plot the frequencies of the number of clusters suggested by
# the internal indexes in the NbClust package
factoextra::fviz_nbclust(nbclust.res, 
   barfill = "#00AFBB", barcolor="#00AFBB", 
   ggtheme = theme_minimal())

# Determining the best initial prototype matrix in ten runs by using 
# the kmpp function for K-means++ algorithm in the inaparc package 
nruns <- 10
miniter <- Inf
for(i in 1:nruns){
  v <- inaparc::kmpp(X, k=nc2)$v
  km.res <- kmeans(X, centers=v, iter.max=100)
  cat(paste0("Prototypes matrix converged ", km.res$iter, " iterations in ", i, "th run.\n"))
  if(km.res$iter < miniter){
    miniter <- km.res$iter
    centers.kmpp <- v
  }
}
cat(paste0("Unsing the best prototypes matrix, built by kmpp initialization,
 K-means converged in ", miniter, " iterations.\n"))
cat("Best prototypes matrix by kmpp is:\n")
print(centers.kmpp)

# K-means clustering with the best initial prototypes matrix
# produced by the kmmp function of the inaparc package
km.res <- kmeans(X, centers=centers.kmpp, nstart = 20)

# Centers of the clusters
cl.medians <- aggregate(X, by=list(cluster=km.res$cluster),median)
cl.means <- aggregate(X, by=list(cluster=km.res$cluster), mean)
cat("Cluster Medians\n")
print(cl.medians)		
cat("Cluster Means\n")
print(cl.means)		

# Or cluster centers from K-means clustering 
cat("Clusters centers from K-means clustering\n")
print(km.res$centers)		

# Visualize the K-means clustering results by the first two features
plot(X, pch=km.res$cluster + 14, col=km.res$cluster)
text((X), col="gray", pos=1, cex=0.7)
points(centers.kmpp, pch=12, col=1:nc2, cex=2.0)
points(km.res$center, pch=10, col=1:nc2, cex=2.0)

# Visualize the K-means clustering results by the feature pairs 
pairs(X, col=km.res$cluster, pch=20)

# Plot the clustering results with clustplot of cluster package
cluster::clusplot(X, km.res$cluster, color=TRUE, 
   shade=TRUE, labels=3, lines=0)

# Aesthetically visualize the clustering results using fviz_cluster of factoextra package
factoextra::fviz_cluster(km.res, data=X,
   ellipse.type = "convex",
   show.clust.cent = TRUE,
   labelsize = 9,
   palette = "jco",
   main = "Cluster plot of K-means clustering with the prototypes initialized by kmpp",
   ggtheme = theme_minimal())

# Generate the initial prototype matrix using the inofrep function of the inaparc package. 
# Since it is based on a deterministic algorithm which produces the same prototype matrix
# in multiple runs, onone run is needed with this prototypes initialization algorithm
centers.inofrep <- inaparc::inofrep(X, k=nc2, nbins=30)$v
inofrep.res <- kmeans(X, centers=centers.inofrep, iter.max=100)

cat("Using the prototypes matrix, built by Inofrep K-means converged in ", inofrep.res$iter, " iterations.\n")
cat("Best prototypes matrix by Inofrep is:\n")
print(centers.inofrep)

# K-means clustering with the best initial prototypes matrix
km.res <- kmeans(X, centers=centers.inofrep, nstart = 20)

# Visualize the K-means clustering results by the first two features
plot(X, pch=km.res$cluster + 14, col=km.res$cluster)
text((X), col="gray", pos=1, cex=0.7)
points(centers.inofrep, pch=12, col=1:nc2, cex=2.0)
points(km.res$center, pch=10, col=1:nc2, cex=2.0)

# Visualize the K-means clustering results by the feature pairs 
pairs(X, col=km.res$cluster, pch=20)

# Plot the clustering results with clustplot of the cluster package
cluster::clusplot(X, km.res$cluster, color=TRUE, 
   shade=TRUE, labels=3, lines=0)

# Aesthetically visualize the clustering results using fviz_cluster of factoextra package
factoextra::fviz_cluster(km.res, data=X,
   ellipse.type = "convex",
   show.clust.cent = TRUE,
   labelsize = 9,
   palette = "jco",
   main = "Cluster plot of K-means clustering with the prototypes initialized by inofrep",
   ggtheme = theme_minimal())

# Visual clustering tendency assessment
gradient.color <- list(low = "#FFFFFF",  high = "#00AFBB")
act <- get_clust_tendency(data=X, n = 10, gradient = gradient.color)
plot(act$plot)
