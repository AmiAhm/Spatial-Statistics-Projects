library(ggplot2)
library(reshape2)
library(latex2exp)
library(geoR)
library(foreach)
library(tidyr)

# Plotting 1a)
from = 1
to = 50
xx <- seq(from = from, to = to, by = 1)
vs.exponential <- c(1.1,1.9)
vs.matern <- c(1, 3)
sigma2s <- c(1, 5)


plot.covariance <- function(cov.model = "matern", kappa = 2, sigma2 = 1, from = 1, to = 50, by = 0.01 , phi = 1, color = NULL){
  x <- seq(from = from ,to = to, by = by)
  y <- cov.spatial(x, cov.pars=c(sigma2, phi), cov.model = cov.model, kappa = kappa)
  p <- geom_line(data = as.data.frame(cbind(x, y)), aes(x = x, y = y, color = color))
  return(p)
}



for_each_paramater_do <- function(fun, ...){
  result = list()
  foreach(sigma2=sigma2s) %:%
    foreach(cov.model = c("matern", "exp")) %do%{
      for(i in 1:2){
        if(cov.model == "matern"){
          vs <- vs.matern[i]
        }else{
          vs <- vs.exponential[i]
        }
        p <- fun(sigma2 = sigma2, cov.model = cov.model, kappa = vs, ...)
        result = c(result, list(p))
      }
    }
    return(result)
}

for_each_parameter_plot <- function(fun, ...){
  plots <- for_each_paramater_do(fun, ...)
  plot.object <- ggplot()
  foreach(p=plots) %do% {plot.object <- plot.object + p}
  return(plot.object)
}
p1 <- for_each_parameter_plot(plot.covariance, from = 1, to = 50, by = 0.01, phi = 10, color = "black")

variogram <- function(x, ...){
  return(1-cov.spatial(x, ...))
}

p2 <- curve(variogram(x, cov.pars=c(1, .2), cov.model = "matern", kappa = 2), from = 0, to = 2,
      xlab = "distance", ylab = expression(gamma(h)),
      main = "Matern Variogram")


# Problem 1b)

# Calculate covariance matrix. Find all combinations of points
library(MASS)

# TODO: Craete some for loop and plot around this
cov.matr <- function(xx, yy, sigma2, phi, cov.model, kappa){
  cov.fun <- function(x) cov.spatial(x, cov.pars=c(sigma2, phi), cov.model = cov.model, kappa = kappa)
  cov.matrix <- matrix(NA, length(xx), length(yy))
  foreach(i=1:length(xx)) %:%
    foreach(j=1:length(yy)) %do%{ # can switch j=1 with j=1 if quadratic
      a <- cov.fun(abs(xx[i]-yy[j]))
      cov.matrix[i, j] = a
    }
  return(cov.matrix)
}

sigma2 <- 1
phi <- 10
cov.model <- "matern"
kappa <- 2
cov.matrix = cov.matr(xx, xx, sigma2, phi, cov.model, kappa)


mu = rep(0, length(xx))
draw <- MASS::mvrnorm(n = 8, mu = mu, Sigma = cov.matrix)
draw <- t(draw)
draw <- cbind(xx, draw)
observations <- Reduce(rbind, lapply(2:ncol(draw), function(col) cbind(draw[,c(1, col)], col)))
observations <- as.data.frame(observations)
colnames(observations) <- c("x", "observed_value", "trial")
observations$trial <- as.factor(observations$trial)

pl <- ggplot()
pl + geom_line(data = observations, aes(x = x, y = observed_value, color = trial))

# Problem 1d) 
sigma2 <- 1
phi <- 10
cov.model <- "matern"
kappa <- 2
dd = c(10, 25, 30)
mu.d = rep(0, length(dd))
mu.l = rep(0, length(xx))
sigma.ll= cov.matr(xx,xx, sigma2, phi, cov.model, kappa)
sigma.dd = cov.matr(dd, dd, sigma2, phi, cov.model, kappa)
sigma.ld = cov.matr(xx, dd, sigma2, phi, cov.model, kappa)
sigma.dl = t(sigma.ld)
sigma.l.d = sigma.ll - sigma.ld %*% solve(sigma.dd) %*% sigma.dl

ddx <- observations[observations$col == 2 & observations$xx %in% dd,]$V2
mu.l.d. = mu.l + sigma.ld %*% solve(sigma.dd) %*% (ddx-mu.d)
# 90 percent confidence:
c = qnorm(p = 0.95)
variances = diag(sigma.l.d) 
variances[variances < 1e-12] = 0
min.90 = mu.l.d. - c*sqrt(variances)
max.90 = mu.l.d. + c*sqrt(variances)
plot(xx, mu.l.d.)
points(xx, min.90)
points(xx, max.90)

