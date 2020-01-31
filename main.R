library(ggplot2)
library(reshape2)
library(latex2exp)
library(geoR)
library(foreach)

# Plotting 1a)
from = 1
to = 50
xx <- seq(from = from, to = to, by = 1)
vs.exponential <- c(1.1,1.9)
vs.matern <- c(1, 3)
sigmas <- c(1, 5)


plot.covariance <- function(cov.model = "matern", kappa = 2, sigma2 = 1, from = 1, to = 50, by = 0.01 , phi = 1){
  x <- seq(from = from ,to = to, by = by)
  y <- cov.spatial(x, cov.pars=c(sigma2, phi), cov.model = cov.model, kappa = kappa)
  p <- ggplot(as.data.frame(cbind(x, y)), aes(x = x, y = y)) + geom_line() + xlab(TeX("$\\tau$")) + ylab(TeX("$\\rho(\\tau)$"))
  return(p)
}

for.parameters.do <- function(param1, param2, param3, fun, ...){
  result <- list()
  foreach(p1 = param1, p2 = param2, p3 = param3)%do% result <- c(result, fun(p1, p2, p3, ...))
  return(result)
}

plot.covariance()

variogram <- function(x, ...){
  return(1-cov.spatial(x, ...))
}

p2 <- curve(variogram(x, cov.pars=c(1, .2), cov.model = "matern", kappa = 2), from = 0, to = 2,
      xlab = "distance", ylab = expression(gamma(h)),
      main = "Matern Variogram")

# Problem 1b)

