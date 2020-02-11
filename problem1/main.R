library(ggplot2)
library(reshape2)
library(latex2exp)
library(geoR)
library(foreach)
library(tidyr)
library("viridis")        
# Plotting 1a)

# Defining some parameters
from = 1
to = 50
xx <- seq(from = from, to = to, by = 1)
vs.exponential <- c(1.1,1.9)
vs.matern <- c(1, 3)
sigma2s <- c(1, 5)

# Creating empty plot frame
colors <- c("blue", "blue", "red", "red","blue", "blue", "red", "red")
x <- seq(from = 0 ,to = 50, by = 0.001)
p1 <- plot(1, type="n", xlab=TeX("$\\tau$"), ylab=TeX("$\\rho_r(\\tau)$"), xlim=c(0, 50), ylim=c(0, 7),
     main = TeX("Matern covariance, $\\sigma^2_r\\rho_r(\\tau)$"))
# Plotting matern stuff
y <- cov.spatial(x, cov.pars=c(sigma2s[1], vs.matern[1]), cov.model = "matern", kappa = 2)
lines(x = x,
      y = y,
      col = colors[1],
      xlim = c(0,50),
      lwd = 2,
      lty = 1)
y <- cov.spatial(x, cov.pars=c(sigma2s[1], vs.matern[2]), cov.model = "matern", kappa = 2)
lines(x = x,
      y = y,
      col = colors[2],
      xlim = c(0,50),
      lwd = 2, 
      lty = 2)
y <- cov.spatial(x, cov.pars=c(sigma2s[2], vs.matern[1]), cov.model = "matern", kappa = 2)
lines(x = x,
      y = y,
      col = colors[3],
      xlim = c(0,50),
      lwd = 2,
      lty = 3)

y <- cov.spatial(x, cov.pars=c(sigma2s[2], vs.matern[2]), cov.model = "matern", kappa = 2)
lines(x = x,
      y = y,
      col = colors[4],
      xlim = c(0,50),
      lwd = 2,
      lty = 4)
legend('topright', legend=TeX(c(paste("$\\nu =$", vs.matern[1], ", $\\sigma^2 =$",sigma2s[1]),
                                paste("$\\nu =$", vs.matern[2], ", $\\sigma^2 =$",sigma2s[1]),
                                paste("$\\nu =$", vs.matern[1], ", $\\sigma^2 =$",sigma2s[2]),
                                paste("$\\nu =$", vs.matern[2], ", $\\sigma^2 =$",sigma2s[2]))),
                                col = colors[1:4],
       lty = c(1,2,3,4),
       lwd = 2, 
       bty="n",
       inset = 0.08)
# Plotting powered exponential
p2 <- plot(1, type="n", xlab=TeX("$\\tau$"), ylab=TeX("$\\sigma^2_r\\rho_r(\\tau)$"), xlim=c(0, 50), ylim=c(0, 7),
     main = TeX("Powered exponential covariance, $\\sigma^2_r\\rho_r(\\tau)$"))
y <- cov.spatial(x, cov.pars=c(sigma2s[1], vs.exponential[1]), cov.model = "exponential", kappa = 2)
lines(x = x,
      y = y,
      col = colors[5],
      xlim = c(0,50),
      lwd = 2,
      lty = 1)
y <- cov.spatial(x, cov.pars=c(sigma2s[1], vs.exponential[2]), cov.model = "exponential", kappa = 2)
lines(x = x,
      y = y,
      col = colors[6],
      xlim = c(0,50),
      lwd = 2, 
      lty = 2)
y <- cov.spatial(x, cov.pars=c(sigma2s[2], vs.exponential[1]), cov.model = "exponential", kappa = 2)
lines(x = x,
      y = y,
      col = colors[7],
      xlim = c(0,50),
      lwd = 2,
      lty =3)

y <- cov.spatial(x, cov.pars=c(sigma2s[2], vs.exponential[2]), cov.model = "exponential", kappa = 2)
lines(x = x,
      y = y,
      col = colors[8],
      xlim = c(0,50),
      lwd = 2,
      lty = 4)
legend('topright', legend=TeX(c(paste("$\\nu =$", vs.matern[1], ", $\\sigma^2 =$",sigma2s[1]),
                                paste("$\\nu =$", vs.matern[2], ", $\\sigma^2 =$",sigma2s[1]),
                                paste("$\\nu =$", vs.matern[1], ", $\\sigma^2 =$",sigma2s[2]),
                                paste("$\\nu =$", vs.matern[2], ", $\\sigma^2 =$",sigma2s[2]))),
       col = colors[5:8],
       lty = c(1,2,3,4),
       lwd = 2,
       inset = 0.08, 
       bty ="n")



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
observation_error <- function(dd, sigma_e_2){
  # dd
  res <- diag(1, length(dd))*sigma_e_2
  res
}

sigma_e_2 <- .25 # Observation error
sigma2 <- 5
phi <- 10
cov.model <- "matern"
kappa <- 2
dd = c(10, 25, 30)
mu.d = rep(0, length(dd))
mu.l = rep(0, length(xx))
sigma.ll= cov.matr(xx,xx, sigma2, phi, cov.model, kappa)
sigma.dd = cov.matr(dd, dd, sigma2, phi, cov.model, kappa) + observation_error(dd, sigma_e_2)
sigma.ld = cov.matr(xx, dd, sigma2, phi, cov.model, kappa)
sigma.dl = t(sigma.ld)
sigma.l.d = sigma.ll - sigma.ld %*% solve(sigma.dd) %*% sigma.dl

# Assume these are what is actually observed
ddy <- observations[observations$trial == 2 & observations$x %in% dd,]$observed_value
mu.l.d. = mu.l + sigma.ld %*% solve(sigma.dd) %*% (ddy-mu.d)
# 90 percent confidence:
c = qnorm(p = 0.95)
variances = diag(sigma.l.d) 
variances[variances < 1e-12] = 0
min.90 = mu.l.d. - c*sqrt(variances)
max.90 = mu.l.d. + c*sqrt(variances)

min.90 <- as.vector(min.90)
max.90 <- as.vector(max.90)


plot(x = xx, y = mu.l.d., ylim = c(-5, 5))
points(x = xx, y = min.90, col = "red")
points(xx, max.90, col = "red")

abline(v=dd[1], lty = 2)
abline(v=dd[2], lty = 2)
abline(v=dd[3], lty = 2)

# Problem 1e) 
sigma_e_2 <- 0 # Observation error
sigma2 <- 5
phi <- 10
cov.model <- "matern"
kappa <- 2
dd = c(10, 25, 30)
mu.d = rep(0, length(dd))
mu.l = rep(0, length(xx))
sigma.ll= cov.matr(xx,xx, sigma2, phi, cov.model, kappa)
sigma.dd = cov.matr(dd, dd, sigma2, phi, cov.model, kappa) + observation_error(dd, sigma_e_2)
sigma.ld = cov.matr(xx, dd, sigma2, phi, cov.model, kappa)
sigma.dl = t(sigma.ld)
sigma.l.d = sigma.ll - sigma.ld %*% solve(sigma.dd) %*% sigma.dl
mu.l.d. = mu.l + sigma.ld %*% solve(sigma.dd) %*% (ddy-mu.d)

c = qnorm(p = 0.95)
variances = diag(sigma.l.d) 
variances[variances < 1e-12] = 0
min.90 = mu.l.d. - c*sqrt(variances)
max.90 = mu.l.d. + c*sqrt(variances)

min.90 <- as.vector(min.90)
max.90 <- as.vector(max.90)

draw <- MASS::mvrnorm(n = 100, mu = mu.l.d., Sigma = sigma.l.d)
draw <- t(draw)
draw <- cbind(xx, draw)
observations <- Reduce(rbind, lapply(2:ncol(draw), function(col) cbind(draw[,c(1, col)], col)))
observations <- as.data.frame(observations)
colnames(observations) <- c("x", "observed_value", "trial")
observations$trial <- as.factor(observations$trial)

plot(x = xx, y = mu.l.d., ylim = c(-5, 5))
points(x = xx, y = min.90, col = "red")
points(xx, max.90, col = "red")
points(observations$x, observations$observed_value, col = rgb(0,1,0,0.25))
abline(v=dd[1], lty = 2)
abline(v=dd[2], lty = 2)
abline(v=dd[3], lty = 2)

## Calculating percentage points inside:
observations$inside_min_max_90 = observations$observed_value > min.90 & observations$observed_value < max.90
mean(observations$inside_min_max_90)

# Problem 1f)

## Chance posterior over 2
plot.height <- 10 
plot.min <- -5
mu <- mu.l.d.
variances <- diag(sigma.l.d)
probs <- pnorm(2, mean = mu, sd = sqrt(variances), lower.tail = F) # TODO: Sqrt here or not?
probs <- probs * plot.height + plot.min


## TODO: Use integrate to find close to theoretical value of A_r 
# Check if point over 2:
observations$over_2 <- observations$observed_value > 2
a_r.est <- mean(observations$over_2)*50 # Points are equally spread out so placement does not matter



# Look at each x
library(dplyr)
grouped_over_2 <- observations %>% group_by(x) %>% summarise(perc_over_2 = mean(over_2))


#### Plotting:
plot(x = xx, y = mu.l.d., ylim = c(-5, 5))
points(x = xx, y = min.90, col = "red")
points(x = xx, y = min.90, col = "red")
points(xx, max.90, col = "red")
points(observations[observations$over_2,]$x, observations[observations$over_2,]$observed_value, col = rgb(0,0,1,0.25))
points(observations[!observations$over_2,]$x, observations[!observations$over_2,]$observed_value, col = rgb(0,1,0,0.25))
abline(v=dd[1], lty = 2)
abline(v=dd[2], lty = 2)
abline(v=dd[3], lty = 2)
lines(xx, probs)
abline(h = 2)
lines(grouped_over_2$x, grouped_over_2$perc_over_2*plot.height + plot.min, col="blue")
sum(grouped_over_2$perc_over_2)
####

# Each sample is a bernoulli trial with p given by x (either blue or black line)
# Each sample thus have variance p(1-p)
# p_hat_i = (Xi1 + ... + Xi100) / 100  
# => var p_hat_i = sum(var(X_i1))/100 = p_hat_i*(1-p_hat_i)
# hat A = sum(p_hat_1, p_hat_2, p_hat_3, ... p_hat_50)
# => var(hat_A) = sum(p_hat_i*(1-p_hat_i))
p.hat.i = grouped_over_2$perc_over_2
var.p.hat.i = p.hat.i * (1-p.hat.i)
var.A.hat = sum(var.p.hat.i)


# Alternative method
A.hat.V2 <- sum(mu.l.d. > 2)

