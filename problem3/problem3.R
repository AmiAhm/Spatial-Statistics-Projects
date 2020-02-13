source("./problem1/main.R")
library(ggplot2)
library(reshape2)
library(latex2exp)
library(geoR)
library(foreach)
library(tidyr)
library(ggpubr)
library(dplyr)
library(fields)
library(RandomFields)
library(akima)
# Problem 3
sigma2 <- 2
xi <- 3

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
Sigma <- create_covariance_matrix(obs.points, obs.points, sigma2 = sigma2, xi = xi)
mu <- rep(0, nrow(obs.points))

# Problem 3a)
# Set seed and draw
set.seed(1)
par(mfrow=c(3,2))

draw1 <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
draw.data1 <- interp(obs.points[,1], obs.points[,2], draw1)
fields::image.plot(draw.data1, main = "a) Realization of discretized RF", xlab = "x", ylab = "y")

fields::image.plot(draw.data1, main = "Realization of discretized RF a) with contours", xlab = "x", ylab = "y")
contour(draw.data1,add=T)


draw2 <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
draw.data2 <- interp(obs.points[,1], obs.points[,2], draw2)
fields::image.plot(draw.data2, main = "b) Realization of discretized RF", xlab = "x", ylab = "y")

fields::image.plot(draw.data2, main = "Realization of discretized RF b) with contours", xlab = "x", ylab = "y")
contour(draw.data2,add=T)

draw3 <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
draw.data3 <- interp(obs.points[,1], obs.points[,2], draw3)
fields::image.plot(draw.data3, main = "b) Realization of discretized RF", xlab = "x", ylab = "y")

fields::image.plot(draw.data3, main = "Realization of discretized RF b) with contours", xlab = "x", ylab = "y")
contour(draw.data3,add=T)

# Problem 3b)
# Finding empirical variogram
library(geoR)
library(ggplot2)
par(mfrow=c(1,1))
M1 <- as.matrix(cbind(obs.points[,1], obs.points[,2], draw1))
geo1 <- as.geodata(M1)
variogram1 <- variog(geo1)
p <- ggplot() + geom_point(aes(x=variogram1$u, y = variogram1$v), colour = "Blue") +
  geom_line(aes(x=variogram1$u, y = variogram1$v), colour = "Blue")

M2 <- as.matrix(cbind(obs.points[,1], obs.points[,2], draw2))
geo2 <- as.geodata(M2)
variogram2 <- variog(geo2)
p <- p + geom_point(aes(x=variogram2$u, y = variogram2$v), colour = "Red")  +
  geom_line(aes(x=variogram2$u, y = variogram2$v), colour = "Red")


M3 <- as.matrix(cbind(obs.points[,1], obs.points[,2], draw3))
geo3 <- as.geodata(M3)
variogram3 <- variog(geo3)
p <- p + geom_point(aes(x=variogram3$u, y = variogram3$v), colour = "Green") +
  geom_line(aes(x=variogram3$u, y = variogram3$v), colour = "Green")


# Finding theoretical variogram
xx <- seq(0, 40, by = 0.1)
yy <- sapply(xx, function(x) (sigma2 - cov(x, 0, xi = xi, sigma2 = 1)))
p <- p + geom_line(aes(x = xx, y = yy))
p <- p + ggtitle("Variograms, theoretical and empirical") + xlab("Distance") + ylab(TeX("$\\gamma_r$"))
p

# Problem 3d)
set.seed(12345)
n = 36
xr <- runif(n)*30
yr <- runif(n)*30
obs <- cbind(xr, yr)
mu <- rep(0, n)
Sigma <- create_covariance_matrix(obs, obs, xi = xi, sigma2 = sigma2)
draw <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
draw.data <- cbind(obs[,1], obs[,2], draw)
geo <- as.geodata(draw.data)
variogram1 <- variog(geo)
p <- ggplot() + geom_point(aes(x=variogram1$u, y = variogram1$v), colour = "Blue") +
  geom_line(aes(x=variogram1$u, y = variogram1$v), colour = "Blue") + ggtitle("Variogram, theoretical and empirical. n=36") + xlab("Distance") + ylab(TeX("$\\gamma_r$"))
p <- p + geom_line(aes(x = xx, y = yy))
p

# Problem 3d)
l.fit <- function(n, p = ggplot(), colour = "Black"){
  set.seed(12345)
  xr <- runif(n)*30
  yr <- runif(n)*30
  obs <- cbind(xr, yr)
  mu <- rep(0, n)
  Sigma <- create_covariance_matrix(obs, obs, xi = xi, sigma2 = sigma2)
  draw <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
  draw.data <- cbind(obs[,1], obs[,2], draw)
  geo <- as.geodata(draw.data)
  # Now estiamte xi and sigma
  mod <- likfit(geo, cov.model = "exponential", ini.cov.pars = c(1,1))
  xi <- mod$phi 
  sigma2 <- mod$sigmasq
  yy <- sapply(xx, function(x) (sigma2 - cov(x, 0, xi = xi, sigma2 = 1)))
  
  variogram1 <- variog(geo)
  p <- p + geom_point(aes(x=variogram1$u, y = variogram1$v), colour = colour) +
    geom_line(aes(x=variogram1$u, y = variogram1$v), colour = colour) +
    geom_line(aes(x = xx, y = yy), colour = colour, linetype = "dashed") +
    geom_line(aes(x = xx, y = yy.theoretical))
  
  return(list(likfit = mod, p = p))
}
yy.theoretical <- sapply(xx, function(x) (sigma2 - cov(x, 0, xi = xi, sigma2 = 1)))

l.fit(36)$likfit
p1 <- l.fit(36, "Blue")$p

# Problem 3e)
p2 <- l.fit(9, p, "Red")$p
p3 <- l.fit(64, p, "Green")$p
p4 <- l.fit(100,p, "Purple")$p
ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2)


