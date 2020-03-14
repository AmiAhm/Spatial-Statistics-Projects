
library(ggplot2)
library(tidyverse)
library(ggforce)
library(spatial)
library(MASS)
library(ggpubr)

# Problem 3
# Plot data
redwood <- read.csv("redwood.dat", header = F, sep = " ")
redwood <- as.data.frame(redwood)
colnames(redwood) <- c("x", "y")

ggplot(data = redwood) + geom_point(aes(x,y)) + theme_classic()


# Doing hierchical clustering with centroid and eculidan
dist_mat <- dist(redwood, method = 'euclidean') # Calculate distances
hclust_centroid <- hclust(dist_mat, method = 'centroid') # Use centroid as moteher node


# Plot data
plot(hclust_centroid, main = "Cluster Dendogram for Redwoods. \n Centroid method w/ euclidian distance")


cut_centroid <- cutree(hclust_centroid, k = 8) # Find node of elements when we have 7 groups



redwood$mother = as.factor(cut_centroid)
ggplot(redwood) + geom_point(aes(x= x, y = y, color = mother, shape = mother)) + theme_classic() +  scale_shape_manual(values=1:nlevels(redwood$mother))


# Finding centers and count for each mother
mothers <- redwood %>% 
  group_by(mother) %>% 
  summarise(x = mean(x),
            y = mean(y), 
            n = n()) %>%
  as.data.frame()


redwood_w_mother <- left_join(redwood, mothers, by = c("mother"))
colnames(redwood_w_mother) <- c("x", "y", "mother", "x.m", "y.m", "n")
redwood_w_mother$dist_mother <- sqrt((redwood_w_mother$x - redwood_w_mother$x.m)^2 + (redwood_w_mother$y - redwood_w_mother$y.m)^2)
redwood_w_mother$x.dist <- redwood_w_mother$x - redwood_w_mother$x.m
redwood_w_mother$y.dist <- redwood_w_mother$y - redwood_w_mother$y.m

std2 <- 1/(length(redwood_w_mother$x.dist) + length(redwood_w_mother$y.dist) - 1)* (sum(redwood_w_mother$x.dist^2) + sum(redwood_w_mother$y.dist^2))

q95 <- qnorm(0.95, mean = 0, sd = sqrt(std2))

ggplot() +
  geom_point(aes(x= redwood$x, y = redwood$y, color = redwood$mother, shape = redwood$mother)) +
  theme_classic() + 
  scale_shape_manual(values=1:nlevels(redwood$mother)) +
  geom_point(aes(x = mothers$x, y = mothers$y)) + 
  geom_circle(aes(r = rep(q95, length(mothers$y)), y0 = mothers$y, x0 = mothers$x), linetype = "dashed", alpha = 0.8)



# Mean of realization around mother
ns <- redwood %>% group_by(mother) %>% summarise(n = n()) %>% ungroup()
lambda.c <- mean(ns$n)



redwood <- read.csv("redwood.dat", header = F, sep = " ")
redwood <- as.data.frame(redwood)
colnames(redwood) <- c("x", "y")

D = ppregion(xl = 0, xu = 1, yl = 0, yu = 1)

redwood.pp2 <- list(
  x = redwood$x, 
  y = redwood$y, 
  area = D
)


kfn <- Kfn(redwood.pp2, 1, k = 200)
plot(kfn, type="b", xlab="distance", ylab="L(t)")

# Distances we want to check out 
neumann_scot_generate <- function(lambda.m, lambda.c, std2, seed){
  set.seed(seed)
  count = 0
  k <- rpois(n = 1, lambda = lambda.m)
  xs <- c()
  mothers <- c()
  for(j in 1:k){
    xj <- runif(2)
    n.child <- rpois(1, lambda = lambda.c)
    for(i in 1:n.child){
      x <- MASS::mvrnorm(n = 1, Sigma = std2*diag(2), mu = xj)
      # Just project for now
      if(x[1] < 0) x[1] <- 0 
      if(x[1] > 1) x[1] <- 1
      
      if(x[2] < 0) x[2] <- 0 
      if(x[2] > 1) x[2] <- 1
      xs <- rbind(xs, x)
    }
    count <- count + n.child
    mothers <- c(mothers, rep(j, n.child))
  }
  list(xs = xs, mothers = mothers) 
}

# Test differnt parameters
lambda.m <- 8
lambda.c <- 7.75
n.rels <- 1000
rels <- lapply(1:n.rels, function(seed) neumann_scot_generate(lambda.m = lambda.m, lambda.c = lambda.c, std2 = std2, seed = seed)$xs)
rels <- lapply(rels, function(rel) list(x = rel[,1], y = rel[,2], area = D))
kfns <- lapply(rels, function(rel) Kfn(rel, 1, k = 200))
bands <- sapply(kfns, function(kfn) kfn$y)
lower_band <- apply(bands, 1, function(x) quantile(x, probs = 0.05))
upper_band <- apply(bands, 1, function(x) quantile(x, probs = 0.95))
mid_band <- apply(bands, 1, function(x) quantile(x, probs = 0.5))

ggplot() + geom_ribbon(aes(x = kfn$x, ymax = upper_band, ymin = lower_band), alpha = 0.3) + geom_line(aes(x = kfn$x, y = kfn$y))

# Seems to be a bit bad a mid distances, a higher variance might on spread might help on that. 
lambda.m <- 7.5
lambda.c <- 6
n.rels <- 1000
rels <- lapply(1:n.rels, function(seed) neumann_scot_generate(lambda.m = lambda.m, lambda.c = lambda.c, std2 = std2, seed = seed)$xs)
rels <- lapply(rels, function(rel) list(x = rel[,1], y = rel[,2], area = D))
kfns <- lapply(rels, function(rel) Kfn(rel, 1, k = 200))
bands <- sapply(kfns, function(kfn) kfn$y)
lower_band <- apply(bands, 1, function(x) quantile(x, probs = 0.05))
upper_band <- apply(bands, 1, function(x) quantile(x, probs = 0.95))
mid_band <- apply(bands, 1, function(x) quantile(x, probs = 0.5))
ggplot() + geom_ribbon(aes(x = kfn$x, ymax = upper_band, ymin = lower_band), alpha = 0.3) + geom_line(aes(x = kfn$x, y = kfn$y))


# Helped somewhat, change lambda somewhat to see if it helps. 
lambda.m <- 7.5
lambda.c <- 6
n.rels <- 1000
std2 <- 1.5*std2
rels <- lapply(1:n.rels, function(seed) neumann_scot_generate(lambda.m = lambda.m, lambda.c = lambda.c, std2 = std2, seed = seed)$xs)
rels <- lapply(rels, function(rel) list(x = rel[,1], y = rel[,2], area = D))
kfns <- lapply(rels, function(rel) Kfn(rel, 1, k = 200))
bands <- sapply(kfns, function(kfn) kfn$y)
lower_band <- apply(bands, 1, function(x) quantile(x, probs = 0.05))
upper_band <- apply(bands, 1, function(x) quantile(x, probs = 0.95))
mid_band <- apply(bands, 1, function(x) quantile(x, probs = 0.5))

ggplot() + geom_ribbon(aes(x = kfn$x, ymax = upper_band, ymin = lower_band), alpha = 0.3) + geom_line(aes(x = kfn$x, y = kfn$y))


# Does not seem to have any good effect.  We try to increase variance again. 
lambda.m <- 7
lambda.c <- 6.8
n.rels <- 1000
std2 <- 2*std2
rels <- lapply(1:n.rels, function(seed) neumann_scot_generate(lambda.m = lambda.m, lambda.c = lambda.c, std2 = std2, seed = seed)$xs)
rels <- lapply(rels, function(rel) list(x = rel[,1], y = rel[,2], area = D))
kfns <- lapply(rels, function(rel) Kfn(rel, 1, k = 200))
bands <- sapply(kfns, function(kfn) kfn$y)
lower_band <- apply(bands, 1, function(x) quantile(x, probs = 0.05))
upper_band <- apply(bands, 1, function(x) quantile(x, probs = 0.95))
mid_band <- apply(bands, 1, function(x) quantile(x, probs = 0.5))

ggplot() + geom_ribbon(aes(x = kfn$x, ymax = upper_band, ymin = lower_band), alpha = 0.3) + geom_line(aes(x = kfn$x, y = kfn$y))


# Display three realiations of the model with this parameters to what was observed (we overlay our grouping to the plots) : 
df1 <- neumann_scot_generate(lambda.m = lambda.m, lambda.c = lambda.c, std2 = std2, seed = 1)
df1 <- as.data.frame(df1)
colnames(df1) <- c("x", "y", "mothers")
df1$mothers <- as.factor(df1$mothers)
df2 <- neumann_scot_generate(lambda.m = lambda.m, lambda.c = lambda.c, std2 = std2, seed = 2)
df2 <- as.data.frame(df2)
colnames(df2) <- c("x", "y", "mothers")
df2$mothers <- as.factor(df2$mothers)
df3 <- neumann_scot_generate(lambda.m = lambda.m, lambda.c = lambda.c, std2 = std2, seed = 3)
df3 <- as.data.frame(df3)
colnames(df3) <- c("x", "y", "mothers")
df3$mothers <- as.factor(df3$mothers)

p1 <- ggplot() +
  geom_point(aes(x= df3$x, y = df3$y, color = df3$mothers, shape = df3$mothers)) +
  theme_classic() + 
  scale_shape_manual(values=1:nlevels(df3$mothers))

p2 <- ggplot() +
  geom_point(aes(x= df2$x, y = df2$y, color = df2$mothers, shape = df2$mothers)) +
  theme_classic() + 
  scale_shape_manual(values=1:nlevels(df2$mothers))

p3 <- ggplot() +
  geom_point(aes(x= df1$x, y = df1$y, color = df1$mothers, shape = df1$mothers)) +
  theme_classic() + 
  scale_shape_manual(values=1:nlevels(df1$mothers))




# Problem 4
# Read and prepare data
cell <- read.csv("cell.dat", header = F, sep = " ")
redwood$mother = as.factor(cut_centroid)
cell <- as.data.frame(cell)
colnames(cell) <- c("x", "y")

ggplot(data = cell) + geom_point(aes(x,y)) + theme_classic()


D = ppregion(xl = 0, xu = 1, yl = 0, yu = 1)

cell.pp <- list(
  x = cell$x, 
  y = cell$y, 
  area = D
)

kfn.cell <- Kfn(cell.pp, 1, k = 200)
plot(kfn.cell)

dist_mat <- dist(cell, method = 'euclidean') # Calculate distances
distances <- as.vector(dist_mat)


ggplot() + geom_density(aes(x=distances), fill = "skyblue", alpha = 0.5) + theme_classic() + xlab("Distance") + ggtitle("Estiamted density of distance between points")

phi <- function(tau, phi0, phi1, tau0){
  if(tau > tau0){
    return(phi0)
  }
  
  phi0*exp(-phi1*abs(tau - tau0))
}

xx <- seq(0, 2, by =0.001)
yy <- sapply(xx, function(x) phi(x, phi0, phi1, tau0))
ggplot() +geom_line(aes(x= xx, y = yy)) + theme_classic() #+ ylim(c(0,1))


tau <- function(x1, x2){
  sqrt(sum((x1-x2)^2)) 
}

gibsstep <- function(df, phi0, tau0, phi1, all = F, ignoreDist = NA){
  u <- runif(nrow(df))
  
  #proposals <- df 
  #proposals$accepted <- FALSE
  #proposals$alpha <- NA
  for(i in 1:nrow(df)){
    xu <- as.vector(df[i,])
    xp <- runif(2)
    
    Res <- apply(df, 1, function(x) {
      left <- tau(xu, x)
      
      left  <- phi(tau = left, phi0 = phi0, phi1 = phi1, tau0 = tau0)
      right <- tau(xp, x)
      if(!is.na(ignoreDist)&& right < ignoreDist){
        return(NA)
      }
      
      right  <- phi(tau = right, phi0 = phi0, phi1 = phi1, tau0 = tau0)
      return(left-right)
    })
    
    res <- -sum(Res)
    alpha <- exp(res)
    alpha <- min(1, alpha)
    #  proposals[i,]$alpha <- alpha
    
    if(!is.na(alpha ) && u[i] < alpha){
      df[i,] = xp
      #   proposals[i,]$accepted <- TRUE
      
    }#else{
    #print(paste(alpha))
    #}
    #proposals[i,c(1,2)] = xp
  }
  
  if(all){
    return(list(df = df, proposals = proposals))
  }
  
  df
}




simulate <- function(n, phi0, tau0, phi1, ignoreDist = NA)
{
  count = 0
  res <- list()
  df.t <- df
  while(count < n){
    print(count)
    df.t <- gibsstep(df.t, phi0, tau0, phi1, ignoreDist = ignoreDist)
    count <- count + 1
    if(count == 1){
      df.1 <- df.t
    }
    #if(count > 20 && count %% 5 == 0)
    res <- c(res, list(df.t))
  }
  res
}

sim <- function(n, phi0, tau0, phi1, ignoreDist = NA){
  res <- simulate(n, phi0 = phi0, tau0 = tau0, phi1 = phi1, ignoreDist = ignoreDist)
  # Creating kfns
  rels <- lapply(res, function(df) list(x = df[,1], y = df[,2], area = D))
  kfns <- lapply(rels, function(rel) Kfn(rel, 1, k = 200))
  bands <- sapply(kfns, function(kfn) kfn$y)
  lower_band <- apply(bands, 1, function(x) quantile(x, probs = 0.05))
  upper_band <- apply(bands, 1, function(x) quantile(x, probs = 0.95))
  mid_band <- apply(bands, 1, function(x) quantile(x, probs = 0.5))
  
  # Distances
  # df <- res[[50]]
  # data <- list(x = df[,1], y = df[,2], area = D)
  # kfn <- Kfn(data, 1, k = 200)
  # plot(kfn)
  # dist_mat <- dist(df, method = 'euclidean') # Calculate distances
  # dist_mat <- as.vector(dist_mat)
  # dist_mat
  
  #Create trace plot
  # 
  # # x and y independent so look at trace for x and ys sepearetly. Look at poistion of first sel
  # pos1.x <- sapply(res, function(df) df[1,1])
  # pos1.y <- sapply(res, function(df) df[1,2])
  # 
  # # Create trace plots 
  # ggplot() + theme_classic() + geom_line(aes(x = 1:length(pos1.x), y = pos1.x))
  # 
  # ggplot() + theme_classic() + geom_line(aes(x = 1:length(pos1.y), y = pos1.y))
  if(is.na(ignoreDist)){
    p <- ggplot() + 
      geom_ribbon(aes(x = kfn.cell$x, ymax = upper_band, ymin = lower_band), alpha = 0.3) +
      geom_line(aes(x = kfn.cell$x, y = kfn.cell$y)) +
      xlab("Distance") + 
      ylab("Kfn") + 
      ggtitle(paste("Simualted cell data \n phi0 =", phi0, ", phi1 =", phi1, ", tau0 =", round(tau0,3))) +
      theme_classic() + theme(text = element_text(size=20))
  }else{
    p <- ggplot() + 
      geom_ribbon(aes(x = kfn.cell$x, ymax = upper_band, ymin = lower_band), alpha = 0.3) +
      geom_line(aes(x = kfn.cell$x, y = kfn.cell$y)) +
      xlab("Distance") + 
      ylab("Kfn") + 
      ggtitle(paste("Simualted cell data \n phi0 =", phi0, ", phi1 =", phi1, ", tau0 =", round(tau0,3), ", tau1 = ", round(ignoreDist,3))) +
      theme_classic() + theme(text = element_text(size=20))
    
  }
  
  return(list(res = res, p  = p))
}

# Simulating Neumm scott for different parameters
set.seed(2)
df <- cell
phi0 <- 1
tau0 <- min(distances)
phi1 <- 1
sim1 <- sim(101, phi0, tau0, phi1)

set.seed(2)
df <- cell
phi0 <- 1
tau0 <- min(distances)
phi1 <- 100
sim2 <- sim(101, phi0, tau0, phi1)



set.seed(2)
df <- cell
phi0 <- 100
tau0 <- min(distances)
phi1 <- 1
sim3 <- sim(101, phi0, tau0, phi1)


set.seed(2)
df <- cell
phi0 <- 50
tau0 <- 0.1
phi1 <- 10000
sim4 <- sim(101, phi0, tau0, phi1)

# Plotting results and saving plots
p <-  ggarrange(sim1$p, sim2$p, sim3$p, sim4$p)
ggsave(filename = "p41.png", p)

# Do not seem to differnt we can plot some realisations
# Plotting some realisations 
p0 <- ggplot(data = df) + geom_point(aes(x,y)) + theme_classic() + ggtitle("Initial realization") + theme(text = element_text(size=20))
p1 <- ggplot(data = sim4$res[[2]]) + geom_point(aes(x,y)) + theme_classic() + ggtitle("First realization") + theme(text = element_text(size=20))
p2 <- ggplot(data = sim4$res[[51]]) + geom_point(aes(x,y)) + theme_classic() + ggtitle("50th realization") + theme(text = element_text(size=20))
p3 <- ggplot(data = sim4$res[[101]]) + geom_point(aes(x,y)) + theme_classic() + ggtitle("100th realization") + theme(text = element_text(size=20))
plot2 <- ggarrange(p0, p1, p2, p3, nrow = 2, ncol = 2)
ggsave(filename = "p42.png", plot2)

# Plot where we ignore dist
set.seed(2)
df <- cell
phi0 <- 1
tau0 <- 0.32
phi1 <- 1
tau1 <- min(distances)
sim1.2 <- sim(101, phi0, tau0, phi1, ignoreDist = tau1)

sim1$p
set.seed(2)
df <- cell
phi0 <- 1
tau0 <- 0.32
phi1 <- 100
tau1 <- min(distances)
sim2.2 <- sim(101, phi0, tau0, phi1, ignoreDist = tau1)



set.seed(2)
df <- cell
phi0 <- 100
tau0 <- 0.32
phi1 <- 1
tau1 <- min(distances)
sim3.2 <- sim(101, phi0, tau0, phi1, ignoreDist = tau1)


set.seed(2)
df <- cell
phi0 <- 50
tau0 <- 0.1
phi1 <- 10000
tau1 <- min(distances)
sim4.2 <- sim(101, phi0, tau0, phi1, ignoreDist = tau1)


# Plotting results and saving plots
p <-  ggarrange(sim1.2$p, sim2.2$p, sim3.2$p, sim4.2$p)
ggsave(filename = "p43.png", p)



p0 <- ggplot(data = df) + geom_point(aes(x,y)) + theme_classic() + ggtitle("Initial realization") + theme(text = element_text(size=20))
p1 <- ggplot(data = sim4.2$res[[2]]) + geom_point(aes(x,y)) + theme_classic() + ggtitle("First realization") + theme(text = element_text(size=20))
p2 <- ggplot(data = sim4.2$res[[51]]) + geom_point(aes(x,y)) + theme_classic() + ggtitle("50th realization") + theme(text = element_text(size=20))
p3 <- ggplot(data = sim4.2$res[[101]]) + geom_point(aes(x,y)) + theme_classic() + ggtitle("100th realization") + theme(text = element_text(size=20))
plot2 <- ggarrange(p0, p1, p2, p3, nrow = 2, ncol = 2)
ggsave(filename = "p44.png", plot2)

