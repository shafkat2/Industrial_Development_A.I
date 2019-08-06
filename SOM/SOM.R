# input Stata file
#library(foreign)
#mydata <- read.dta("bangla.dta")

#library(haven)
#data <- read_dta('bangla.dta')

#library(readstata13)
#tempdata <- read.dta13('bangla.dta')

#library(readxl)
#library(tseries)
my_data <- read.csv("GyroLinearAcceleration-2019-1-6-6-12-38.csv")
#arma(my_data, order = c(1, 1), lag = NULL, coef = NULL, include.intercept = TRUE, series = NULL, qr.tol = 1e-07)

# Load the kohonen package 
require(kohonen)

# Create a training data set (rows are samples, columns are variables
# Here I am selecting a subset of my variables available in "data"
data_train <- my_data[, c(3,4,5,6,7,8)]

# Change the data frame with training data to a matrix
# Also center and scale all variables to give them equal importance during
# the SOM training process. 
data_train_matrix <- as.matrix(scale(data_train))

# Create the SOM Grid - you generally have to specify the size of the 
# training grid prior to training the SOM. Hexagonal and Circular 
# topologies are possible
som_grid <- somgrid(xdim = 73, ydim=73, topo="hexagonal")

# Finally, train the SOM, options for the number of iterations,
# the learning rates, and the neighbourhood are available
som_model <- som(data_train_matrix, grid=som_grid, rlen=100, alpha=c(0.05,0.01), keep.data = TRUE)


plot(som_model, type="changes")
plot(som_model, type="count")
plot(som_model, type="dist.neighbours")
plot(som_model, type="codes")
coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}
plot(som_model, type = "property", property = som_model$codes[[1]], main=names(som_model$data)[1:5329,1:6], palette.name=coolBlueHotRed)


var <- 6 #define the variable to plot 
var_unscaled <- aggregate(as.numeric(data_train[,var]), by=list(som_model$unit.classif), FUN=mean, simplify=TRUE)[,6] 
plot(som_model, type = "property", property=var_unscaled, main=names(data_train)[var], palette.name=coolBlueHotRed)



mydata <- som_model$codes[[1]] 
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var)) 
for (i in 2:15) {
  wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
}
plot(wss)

## use hierarchical clustering to cluster the codebook vectors
som_cluster <- cutree(hclust(dist(som_model$codes[[1]])), 6)
# plot these results:
pretty_palette <- c("#1f77b4","#ff7f0e","#2ca02c", "#d62728","#9467bd","#8c564b","#e377c2")
plot(som_model, type="mapping", bgcol = pretty_palette[som_cluster], main = "Clusters") 
add.cluster.boundaries(som_model, som_cluster)


# newdata <- read.csv("RMSE.csv")
