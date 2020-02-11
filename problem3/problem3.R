source("./problem1/main.R")

# Problem 1
obs.points <- expand.grid(x = 1:30, y =1:30)
obs.points <- as.matrix(grid)
dist <- function(x, y){
  sqrt(sum((x-y)^2))
}

cor <- function(x,y, xi, sigma2){
  tau <- dist(x, y)
  sigma2 * exp(-tau/xi)
}


create_covariance_matrix <- function(X, Y, xi, sigma2){
  res <- matrix(NA, nrow=nrow(X), ncol = nrow(Y))
  for(i in 1:nrow(X)){
    for(j in 1:nrow(Y)){
      res[i, j] = cor(X[i,], Y[j,], xi, sigma2)
    }
  }
  res
}

Sigma <- create_covariance_matrix(obs.points, obs.points, 2, 3)
mu <- rep(0, nrow(obs.points))


draw <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)

library(fields)
fields::image.plot(x = obs.points[,1], y = obs.points[,2], z= draw)

x <- 1:30
y <- 1:30

obs.grid <- matrix(NA, nrow = nrow(obs.points),ncol = nrow(obs.points))

length.y <- 30
length.x <- 30
for(i in 1:length.x){
  for(j in 1:length.y){
    i <- as.numeric(i)
    j <- as.numeric(j)
    print(((j-1)*length.y)+i)
    obs.grid[i,j] <- draw[((j-1)*length.y)+i]
  }
}
