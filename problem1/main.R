library(ggplot2)
library(reshape2)
library(latex2exp)
library(geoR)
library(foreach)
library(tidyr)
library(ggpubr)
library(dplyr)
library(RandomFields)
library(MASS)
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


n.realizations <- function(n, xx, sigma2, phi, cov.model, kappa = 2, title = "title"){
  set.seed(1)
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
sigma2 <- 5
phi <- vs.matern[2]
cov.model <- "matern"
kappa <- 2

# Structure observations from 
set.seed(1234)
cov.matrix = cov.matr(xx, xx,sigma2, phi, cov.model = cov.model, kappa)
mu = rep(0, length(xx))
# Draw from multivariate normal
draw <- MASS::mvrnorm(n = n, mu = mu, Sigma = cov.matrix)
observations <-draw[1,]

# Function that retruns ovservations error id matrix
observation_error <- function(dd, sigma_e_2){
  # dd
  res <- diag(1, length(dd))*sigma_e_2
  res
}
# Observations points
dd = c(10, 25, 30)

# Assume these are what is actually observed
ddy <-observations[dd]


generate_plot <- function(title, sigma_e_2){
  # Generate conditional paramteres
  mu.d = rep(0, length(dd))
  mu.r = rep(0, length(xx))
  sigma.rr= cov.matr(xx,xx, sigma2, phi, cov.model, kappa)
  
  sigma.dd = cov.matr(dd, dd, sigma2, phi, cov.model, kappa) + observation_error(dd, sigma_e_2)
  sigma.rd = cov.matr(xx, dd, sigma2, phi, cov.model, kappa)
  sigma.dr = t(sigma.rd)
  sigma.r.d = sigma.rr - sigma.rd %*% solve(sigma.dd) %*% sigma.dr
  mu.r.d. = mu.r + sigma.rd %*% solve(sigma.dd) %*% (ddy-mu.d)
  # 90 percent confidence:
  c = qnorm(p = 0.95)
  variances = diag(sigma.r.d) 
  variances[variances < 1e-12] = 0
  min.90 = mu.r.d. - c*sqrt(variances)
  max.90 = mu.r.d. + c*sqrt(variances)
  
  min.90 <- as.vector(min.90)
  max.90 <- as.vector(max.90)
  
  p <- ggplot(data = as.data.frame(cbind(xx, mu.r.d., min.90, max.90, observations))) +
    geom_line(aes(x = xx, y = mu.r.d.)) +
    geom_line(aes(x = xx, y = min.90), linetype = "dashed") + 
    geom_line(aes(x = xx, y =max.90), linetype = "dashed") +
    geom_ribbon(aes(x = xx, ymin = min.90, ymax = max.90), alpha = 0.5, fill = "skyblue") + 
    geom_point(aes(x = xx, y = mu.r.d.)) +
    geom_point(aes(x = xx, y = min.90), alpha = 0.5) +
    geom_point(aes(x = xx, y = max.90), alpha = 0.5) +
    geom_point(aes(x = xx, y = observations), alpha = 0.5, color = "red") +
    geom_line(aes(x = xx, y =observations), linetype = "dashed", color = "red") +
    geom_vline(xintercept = dd) + 
    ggtitle(title) +
    xlab("x") +
    ylab("r|d") 
    
  p
}
p1 <- generate_plot("a) 0.9 confidence without observation error", 0)
p2 <- generate_plot("b) 0.9 confidence with observation error", .25)
  
ggarrange(p1,p2, ncol=2, widths = c(1,1), heights = c(1,1), nrow = 1)

  


# Problem 1e) 
generate_plot2 <- function(title, sigma_e_2){
  mu.d = rep(0, length(dd))
  mu.r = rep(0, length(xx))
  sigma.rr= cov.matr(xx,xx, sigma2, phi, cov.model, kappa)
  
  sigma.dd = cov.matr(dd, dd, sigma2, phi, cov.model, kappa) + observation_error(dd, sigma_e_2)
  sigma.rd = cov.matr(xx, dd, sigma2, phi, cov.model, kappa)
  sigma.dr = t(sigma.rd)
  sigma.r.d = sigma.rr - sigma.rd %*% solve(sigma.dd) %*% sigma.dr
  mu.r.d. = mu.r + sigma.rd %*% solve(sigma.dd) %*% (ddy-mu.d)
  
  # 90 percent confidence:
  c = qnorm(p = 0.95)
  variances = diag(sigma.r.d) 
  variances[variances < 1e-12] = 0
  min.90 = mu.r.d. - c*sqrt(variances)
  max.90 = mu.r.d. + c*sqrt(variances)
  
  min.90 <- as.vector(min.90)
  max.90 <- as.vector(max.90)
  set.seed(1234)
  draw <- MASS::mvrnorm(n = 100, mu = mu.r.d., Sigma = sigma.r.d)
  draw <- t(draw)
  observations <- cbind(xx, observations)
  observations <- Reduce(rbind, lapply(2:ncol(observations), function(col) cbind(draw[,c(1, col)], col)))
  observations <- as.data.frame(observations)
  colnames(observations) <- c("x", "observed_value", "trial")
  observations$trial <- as.factor(observations$trial)
  
  # Calculate empirical bands
  observations <- observations %>% group_by(x) %>% summarize(emp_mu = mean(observed_value), emp.min_90 = quantile(observed_value, probs = c(0.05)), emp.max_90 = quantile(observed_value, probs = c(0.95)))
  
  # Add theoretical 090% bands. 
  observations$min.90 <- min.90
  observations$max.90 <- max.90 
  observations$mu.r.d. <- mu.r.d. 
  
  # Create plot
  p <- ggplot(data = observations) +
    
    # Empiricla values
    geom_line(aes(x = x, y = emp_mu)) +
    geom_line(aes(x = x, y = emp.min_90), linetype = "dashed") + 
    geom_line(aes(x = x, y =emp.max_90), linetype = "dashed") +
    geom_ribbon(aes(x = x, ymin = emp.min_90, ymax = emp.max_90), alpha = 0.5, fill = "skyblue") + 
    geom_point(aes(x = x, y = emp_mu)) +
    geom_point(aes(x = x, y = emp.min_90), alpha = 0.5) +
    geom_point(aes(x = x, y = emp.max_90), alpha = 0.5) +
    
    # Theoretical values
    geom_line(aes(x = xx, y = mu.r.d.)) +
    geom_line(aes(x = xx, y = min.90), linetype = "dashed") + 
    geom_line(aes(x = xx, y =max.90), linetype = "dashed") +
    geom_ribbon(aes(x = xx, ymin = min.90, ymax = max.90), alpha = 0.5, fill = "red") + 
    geom_point(aes(x = xx, y = mu.r.d.)) +
    geom_point(aes(x = xx, y = min.90), alpha = 0.5) +
    geom_point(aes(x = xx, y = max.90), alpha = 0.5) +
    
    geom_vline(xintercept = dd) + 
    xlab("x") +
    ylab("r|d") +
    ggtitle(title)
  
  list(p=p, draw = draw, mu.r.d.= mu.r.d.)
}
p1 <- generate_plot2("a) 100 realizations, emprical estimation. No observation error", 0)$p
p2 <- generate_plot2("b) 100 realizations, emprical estimation. With observation error", 0.25)$p
ggarrange(p1,p2, ncol=2, widths = c(1,1), heights = c(1,1), nrow = 1)

# Problem 1f)
# Estiamte that each realisation line follow straight line of observed
## Chance posterior over 2
draw <- generate_plot2("", 0)$draw
est_int <- apply(draw, 2, function(x) length(x[x>2]))
est_int
mean(est_int)
var(est_int)

mu.d = rep(0, length(dd))
mu.r = rep(0, length(xx))
sigma.rr= cov.matr(xx,xx, sigma2, phi, cov.model, kappa)

sigma.dd = cov.matr(dd, dd, sigma2, phi, cov.model, kappa)
sigma.rd = cov.matr(xx, dd, sigma2, phi, cov.model, kappa)
sigma.dr = t(sigma.rd)
sigma.r.d = sigma.rr - sigma.rd %*% solve(sigma.dd) %*% sigma.dr
mu.r.d. = mu.r + sigma.rd %*% solve(sigma.dd) %*% (ddy-mu.d)

# 90 percent confidence:
c = qnorm(p = 0.95)
variances = diag(sigma.r.d) 
variances[variances < 1e-12] = 0
min.90 = mu.r.d. - c*sqrt(variances)
max.90 = mu.r.d. + c*sqrt(variances)

min.90 <- as.vector(min.90)
max.90 <- as.vector(max.90)

plot.height <- 10 
plot.min <- -5
mu <- mu.r.d.
variances <- diag(sigma.r.d)
probs <- pnorm(2, mean = mu, sd = sqrt(variances), lower.tail = F) # TODO: Sqrt here or not?
probs <- probs * plot.height + plot.min

# Look at each x
library(dplyr)

observations <- cbind(1:50, draw)
observations <- lapply(1:100, function(i) cbind(observations[,1], observations[,i+1], i))
observations <- Reduce(rbind, observations)
observations <- as.data.frame(observations)
observations$over_2 <- observations$V2>2
grouped_over_2 <- observations %>% group_by(V1) %>% summarise(perc_over_2 = mean(over_2))


#### Plotting:
plot(x = xx, y = mu.r.d., ylim = c(-5, 5))
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
lines(grouped_over_2$V1, grouped_over_2$perc_over_2*plot.height + plot.min, col="blue")
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

require(pracma)
probs <- pnorm(2, mean = mu, sd = sqrt(variances), lower.tail = F) # TODO: Sqrt here or not?
trapz(xx, probs)
trapz(xx, grouped_over_2$perc_over_2)


# Alternative method
A.hat.V2 <- sum(mu.l.d. > 2)

