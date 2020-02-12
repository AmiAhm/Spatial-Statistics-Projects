library(ggplot2)
library(reshape2)
library(latex2exp)
library(geoR)
library(foreach)
library(tidyr)
library(ggpubr)
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
#colors <- c("blue", "blue", "red", "red","blue", "blue", "red", "red")
x <- seq(from = 0 ,to = 50, by = 0.001)
# Plotting matern stuff
y <- cov.spatial(x, cov.pars=c(sigma2s[1], vs.matern[1]), cov.model = "matern", kappa = 2)
p1 <- ggplot() + geom_line(aes(x = x, y= y,col = "red"), data = as.data.frame(cbind(x,y)))
y <- cov.spatial(x, cov.pars=c(sigma2s[1], vs.matern[2]), cov.model = "matern", kappa = 2)
p1 <- p1 + geom_line(aes(x = x, y= y, col = "blue"), data = as.data.frame(cbind(x,y)))
p1 <- p1 + ggtitle(TeX("Matern correlation, $\\rho_r(\\tau)$")) + xlab(TeX("$\\tau$")) + ylab(TeX("$\\rho_r(\\tau)$"))
p1 <- p1 +
  scale_color_discrete("",labels = unname(TeX(c(paste("$\\nu =$", vs.matern[1])
                                             ,paste("$\\nu =$", vs.matern[2]))))) +
  theme_classic() +
  theme(legend.key.size = unit(1.5, 'lines'),
        legend.position = 'right',
        text = element_text(size=20),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"))
  
p1

# Plotting powered exponential
y <- cov.spatial(x, cov.pars=c(sigma2s[1], vs.exponential[1]), cov.model = "exponential", kappa = 2)
p2 <- ggplot() + geom_line(aes(x = x, y= y,col = "red"), data = as.data.frame(cbind(x,y)))
y <- cov.spatial(x, cov.pars=c(sigma2s[1], vs.exponential[2]), cov.model = "matern", kappa = 2)
p2 <- p2 + geom_line(aes(x = x, y= y, col = "blue"), data = as.data.frame(cbind(x,y)))
p2 <- p2 + ggtitle(TeX("Powered exponential correlation, $\\rho_r(\\tau)$")) + xlab(TeX("$\\tau$")) + ylab(TeX("$\\rho_r(\\tau)$"))
p2 <- p2 +
  scale_color_discrete("",labels = unname(TeX(c(paste("$\\nu =$", vs.exponential[1])
                                                ,paste("$\\nu =$", vs.exponential[2])))),
                       guide = guide_legend(label.hjust = 0.1)) +
  theme_classic() +
  theme(legend.key.size = unit(1.5, 'lines'),
        legend.position = 'right',
        text = element_text(size=20),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"))
p2

# variogram <- function(x, ...){
#   return(1-cov.spatial(x, ...))
# }

### Variograms
y <- sigma2s[1]*(1- cov.spatial(x, cov.pars=c(1, vs.matern[1]), cov.model = "matern", kappa = 2))
p3 <- ggplot() + geom_line(aes(x = x, y= y,col = "red"), data = as.data.frame(cbind(x,y)))
y <- sigma2s[1]*(1- cov.spatial(x, cov.pars=c(1, vs.matern[2]), cov.model = "matern", kappa = 2))
p3 <- p3 + geom_line(aes(x = x, y= y,col = "red"), data = as.data.frame(cbind(x,y)), linetype = "dashed")
y <- sigma2s[2]*(1- cov.spatial(x, cov.pars=c(1, vs.matern[1]), cov.model = "matern", kappa = 2))
p3 <- p3 + geom_line(aes(x = x, y= y,col = "blue"), data = as.data.frame(cbind(x,y)))
y <- sigma2s[2]*(1- cov.spatial(x, cov.pars=c(1, vs.matern[2]), cov.model = "matern", kappa = 2))
p3 <- p3 + geom_line(aes(x = x, y= y,col = "blue"), data = as.data.frame(cbind(x,y)), linetype = "dashed")
p3 <- p3 + ggtitle(TeX("Matern variogram, $\\gamma_r(\\tau)$")) + xlab(TeX("$\\tau$")) + ylab(TeX("$\\gamma_r(\\tau)$")) +
  scale_color_discrete("",labels = unname(TeX(c(paste("$\\nu =$", vs.matern[1], ", $\\sigma^2 =$",sigma2s[1])
                                                ,paste("$\\nu =$", vs.matern[2], ", $\\sigma^2 =$",sigma2s[1])
                                                ,paste("$\\nu =$", vs.matern[1], ", $\\sigma^2 =$",sigma2s[2])
                                                ,paste("$\\nu =$", vs.matern[2], ", $\\sigma^2 =$",sigma2s[2])))),
                       guide = guide_legend(label.hjust = 0.1)) +
  theme_classic() +
  theme(legend.key.size = unit(1.5, 'lines'),
        legend.position = 'right',
        text = element_text(size=20),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"))
p3


# Powered exponential
y <- sigma2s[1]*(1- cov.spatial(x, cov.pars=c(1, vs.exponential[1]), cov.model = "exponential", kappa = 2))
p4 <- ggplot() + geom_line(aes(x = x, y= y,col = "red"), data = as.data.frame(cbind(x,y)))
y <- sigma2s[1]*(1- cov.spatial(x, cov.pars=c(1, vs.exponential[2]), cov.model = "exponential", kappa = 2))
p4 <- p4 + geom_line(aes(x = x, y= y,col = "red"), data = as.data.frame(cbind(x,y)), linetype = "dashed")
y <- sigma2s[2]*(1- cov.spatial(x, cov.pars=c(1, vs.exponential[1]), cov.model = "exponential", kappa = 2))
p4 <- p4 + geom_line(aes(x = x, y= y,col = "blue"), data = as.data.frame(cbind(x,y)))
y <- sigma2s[2]*(1- cov.spatial(x, cov.pars=c(1, vs.exponential[2]), cov.model = "exponential", kappa = 2))
p4 <- p4 + geom_line(aes(x = x, y= y,col = "blue"), data = as.data.frame(cbind(x,y)), linetype = "dashed")
p4 <- p4 + ggtitle(TeX("Powered exponential variogram, $\\gamma_r(\\tau)$")) + xlab(TeX("$\\tau$")) + ylab(TeX("$\\gamma_r(\\tau)$")) +
  scale_color_discrete("",labels = unname(TeX(c(paste("$\\nu =$", vs.exponential[1], ", $\\sigma^2 =$",sigma2s[1])
                                                ,paste("$\\nu =$", vs.exponential[2], ", $\\sigma^2 =$",sigma2s[1])
                                                ,paste("$\\nu =$", vs.exponential[1], ", $\\sigma^2 =$",sigma2s[2])
                                                ,paste("$\\nu =$", vs.exponential[2], ", $\\sigma^2 =$",sigma2s[2])))),
                       guide = guide_legend(label.hjust = 0.1)) +
  theme_classic() +
  theme(legend.key.size = unit(1.5, 'lines'),
        legend.position = 'right',
        text = element_text(size=20),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"))
p4




# Problem 1b)

# Calculate covariance matrix. Find all combinations of points
library(MASS)

# Function that returns covariance matrix using cov.spatial
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

n.realizations <- function(n, xx, sigma2, phi, cov.model, kappa = 2, title = "title"){
  # Create covariance matrix
  cov.matrix = cov.matr(xx, xx, sigma2, phi, cov.model, kappa)
  
  # Expected value
  mu = rep(0, length(xx))
  
  # Draw from multivariate normal
  draw <- MASS::mvrnorm(n = n, mu = mu, Sigma = cov.matrix)
  
  # Transform the data
  draw <- t(draw)
  draw <- cbind(xx, draw)
  observations <- Reduce(rbind, lapply(2:ncol(draw), function(col) cbind(draw[,c(1, col)], col)))
  observations <- as.data.frame(observations)
  colnames(observations) <- c("x", "observed_value", "trial")
  observations$trial <- as.factor(observations$trial)
  
  # Plot the data
  pl <- ggplot()
  pl <- pl + geom_line(data = observations, aes(x = x, y = observed_value, color = trial)) + ggtitle(TeX(paste(title, "$\\nu =$", phi, ", $\\sigma^2 =$",sigma2)))
  pl
}

n <- 4
p1 <- n.realizations(n, xx, sigma2s[1], vs.matern[1], cov.model = "matern", kappa = 2, title = "a) Matern realizations, with:")
p2 <- n.realizations(n, xx, sigma2s[1], vs.matern[2], cov.model = "matern", kappa = 2, title = "b) Matern realizations, with:")
p3 <- n.realizations(n, xx, sigma2s[2], vs.matern[1], cov.model = "matern", kappa = 2, title = "c) Matern realizations, with:")
p4 <- n.realizations(n, xx, sigma2s[2], vs.matern[2], cov.model = "matern", kappa = 2, title = "d) Matern realizations, with:")
p5 <- n.realizations(n, xx, sigma2s[1], vs.exponential[1], cov.model = "exponential", kappa = 2, title = "e) Powered exponential realizations, with:")
p6 <- n.realizations(n, xx, sigma2s[1], vs.exponential[2], cov.model = "exponential", kappa = 2, title = "f) Powered exponential, with:")
p7 <- n.realizations(n, xx, sigma2s[2], vs.exponential[1], cov.model = "exponential", kappa = 2, title = "g) Powered exponential, with:")
p8 <- n.realizations(n, xx, sigma2s[2], vs.exponential[2], cov.model = "exponential", kappa = 2, title = "h) Powered exponential, with:")
ggarrange(p1,p2, p3, p4, p5, p6, p7, p8, ncol=2, widths = c(1,1), heights = c(1,1), nrow = 4)



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

