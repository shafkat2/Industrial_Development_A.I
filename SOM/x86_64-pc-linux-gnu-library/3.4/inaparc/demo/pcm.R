# Demo for Possibilistic c-means clustering using kmpp and inofrep 
# prototypes initialization functions in the package 'inaparc'
#

# Install and/or load the required libraries
for(package in c("inaparc", "kpeaks", "cluster", "factoextra", "NbClust", "vegclust")){
  if(!require(package, character.only=TRUE, quietly=TRUE)){
    install.packages(package)
  }
  library(package, character.only=TRUE)
}

# Function to compute the reference distance vector (eta), is used by PCM
cveta <- function(x, t, K=1, m=2){
  veta <- c()
  for(i in 1:ncol(t)){
    total1 <- 0 ; total2 <- 0
    for(j in 1:nrow(t)){
      total1 <- total1 + (t[j,i]^m) * (sum((x[j, ] - v[i, ])^2))
      total2 <- total2 + t[j,i]^m
    }
    veta[i] <- K * total1 / total2 
  }
  return(veta)
}

# Data loading and pre-processing

# Use Iris data set
data(iris)

# Pick the features except the 5th column 
X <- iris[,-5]

# Convert type of X to the matrix class
X <- as.matrix(X)
str(X)
head(X)
tail(X)

# Omit the NAs, if any
X <- na.omit(X)

# For standardization, uncomment the following two lines
# X <- scale(X, center = TRUE, scale = TRUE) 
# head(X); tail(X)

par(ask=TRUE)
# Compute the distance matrix of X
dist.x <- factoextra::get_dist(X, stand = TRUE, method = "pearson")

# Visualize the distance matrix 
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
  cat("K-means converged ", km.res$iter, " iterations in ", i, "th run.\n")
  if(km.res$iter < miniter){
    miniter <- km.res$iter
    centers.kmpp <- v
  }
}
cat("Best prototypes matrix, built by kmpp initialization, K-means converged in ", miniter, " iterations.\n")
cat("Best prototypes matrix by kmpp is:\n")
print(centers.kmpp)

# Firstly run FCM with the best initial prototypes matrix
# created by the kmmp initialization function of the inaparc package
fcm.res <- vegclust::vegclust(X, method="FCM", 
   mobileCenters=centers.kmpp,
   m=2, iter.max=100)

# Compute eta vector to use by PCM
veta <- cveta(X, fcm.res$dist2clusters)

# Use the membership degrees from the previous run of Fuzzy c-means clustering 
centers.kmpp1 <- fcm.res$mobileCenters
pcm.res <- vegclust::vegclust(X, 
   method="PCM", mobileCenters=centers.kmpp1, 
   m=2, eta=veta, iter.max=100)

# Assign hard cluster labels to the objects
clusters <- c(); memberships <- as.matrix(fcm.res$memb)
for(i in 1:nrow(memberships))
  clusters[i] <- which(max(memberships[i,])==memberships[i,])

# Visualize the PCM clustering results by the first two features
plot(X, pch=clusters + 14, col=clusters)
text((X), col="gray", pos=1, cex=0.7)
points(centers.kmpp1, pch=12, col=1:nc2, cex=2.0)
points(pcm.res$mobileCenters, pch=10, col=1:nc2, cex=2.0)

# Visualize the K-means clustering results by the feature pairs 
pairs(X, col=clusters, pch=20)

# Plot the clustering results with clustplot of cluster package
cluster::clusplot(X, clusters, color=TRUE, shade=TRUE, labels=3, lines=0)

# Generate the initial prototype matrix using the inofrep function in the inaparc package. 
# Since it is based on a deterministic algorithm which produces the same prototype matrix
# in multiple runs, onone run is needed with this prototypes initialization algorithm
centers.inofrep <- inaparc::inofrep(X, k=nc2, nbins=30)$v
inofrep.res <- kmeans(X, centers=centers.inofrep, iter.max=100)

cat("Using the prototypes matrix, built by Inofrep, K-means converged in ", inofrep.res$iter, " iterations.\n")
cat("Best prototypes matrix by Inofrep is:\n")
print(centers.inofrep)

# Firstly, run FCM with the best initial prototypes matrix
fcm.res <- vegclust::vegclust(X, method="FCM", 
    mobileCenters=centers.inofrep, 
    m=2, iter.max=100)

# Compute eta vector to use in PCM
veta <- cveta(X, fcm.res$dist2clusters)

# Run PCM with the membership degrees returned from the previous run of FCM
centers.inofrep1 <- fcm.res$mobileCenters
pcm.res <- vegclust::vegclust(X, method="PCM", 
    mobileCenters=centers.inofrep1, 
    m=2, eta=veta, iter.max=100)

# Assign hard cluster labels to the objects
clusters <- c(); memberships <- as.matrix(fcm.res$memb)
for(i in 1:nrow(memberships))
  clusters[i] <- which(max(memberships[i,])==memberships[i,])

# Visualize the PCM clustering results by the first two features 
plot(X, pch=clusters + 14, col=clusters)
text((X), col="gray", pos=1, cex=0.7)
points(centers.inofrep, pch=12, col=1:nc2, cex=2.0)
points(fcm.res$mobileCenters, pch=10, col=1:nc2, cex=2.0)

# Visualize the K-means clustering results by the feature pairs 
pairs(X, col=clusters, pch=20)

# Plot the clustering results with clustplot of cluster package
cluster::clusplot(X, clusters, color=TRUE, 
   shade=TRUE, labels=3, lines=0)

# Visual clustering tendency assessment
gradient.color <- list(low = "#FFFFFF",  high = "#00AFBB")
act <- get_clust_tendency(data=X, n = 10, gradient = gradient.color)
plot(act$plot)
