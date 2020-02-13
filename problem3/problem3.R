source("./problem1/main.R")

# Problem 3

# Craeting grid
obs.points <- expand.grid(x = 1:30, y =1:30)
obs.points <- as.matrix(obs.points)
dist <- function(x, y){
  sqrt(sum((x-y)^2))
}

# Cov function
cov <- function(x,y, xi, sigma2){
  tau <- dist(x, y)
  sigma2 * exp(-tau/xi)
}


# Create covariance matrix
create_covariance_matrix <- function(X, Y, xi, sigma2){
  res <- matrix(NA, nrow=nrow(X), ncol = nrow(Y))
  for(i in 1:nrow(X)){
    for(j in 1:nrow(Y)){
      res[i, j] = cov(X[i,], Y[j,], xi, sigma2)
    }
  }
  res
}

# Create covariance matrinx, and ,u
Sigma <- create_covariance_matrix(obs.points, obs.points, 2, 3)
mu <- rep(0, nrow(obs.points))

# Problem 3a)
# Set seed and draw
set.seed(1)
par(mfrow=c(3,2))

draw1 <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
draw.data1 <- interp(obs.points[,1], obs.points[,2], draw1)
fields::image.plot(draw.data1, main = "a) Realization of discretized RF", xlab = "x", ylab = "y")

fields::image.plot(draw.data1, main = "Realization of discretized RF a) with contours", xlab = "x", ylab = "y")
contour(draw.data,add=T)


draw2 <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
draw.data2 <- interp(obs.points[,1], obs.points[,2], draw2)
fields::image.plot(draw.data, main = "b) Realization of discretized RF", xlab = "x", ylab = "y")

fields::image.plot(draw.data2, main = "Realization of discretized RF b) with contours", xlab = "x", ylab = "y")
contour(draw.data2,add=T)

draw3 <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
draw.data3 <- interp(obs.points[,1], obs.points[,2], draw3)
fields::image.plot(draw.data3, main = "b) Realization of discretized RF", xlab = "x", ylab = "y")

fields::image.plot(draw.data3, main = "Realization of discretized RF b) with contours", xlab = "x", ylab = "y")
contour(draw.data,add=T)

# Problem 3b)
# Finding empirical variogram
library(geoR)
par(mfrow=c(1,1))
M1 <- as.matrix(cbind(obs.points[,1], obs.points[,2], draw1))
geo1 <- as.geodata(M1)
variogram1 <- variog(geo1)
plot(x=variogram1$u, y = variogram$v)

M2 <- as.matrix(cbind(obs.points[,1], obs.points[,2], draw2))
geo2 <- as.geodata(M2)
variogram2 <- variog(geo2)
plot(variogram2, add = T)


# Finding theoretical variogram
sigma2 <- 2
xi <- 3
xx <- seq(0, 40, by = 0.1)
yy <- sapply(xx, function(x) (sigma2 - cov(x, 0, xi, 1)))
points(x =xx, y = yy)

ggplot() + ggim 


# Problem 3c)
n = 36
xr <- runif(n)*30
yr <- runif(n)*30
obs <- cbind(xr, yr)
mu <- rep(0, n)
Sigma <- create_covariance_matrix(obs, obs, xi = xi, sigma2 = sigma2)
draw <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
draw.data <- cbind(obs[,1], obs[,2], draw)
geo <- as.geodata(draw.data)
v1 <- variog(geo)
plot(v1)
