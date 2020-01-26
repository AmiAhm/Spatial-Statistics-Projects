library(ggplot2)
library(reshape2)
library(latex2exp)

tau_10 <- function(X, Y){
  return(abs(X-Y)/10)
}

get_powered_exponential <- function(v){
  if(v < 1 | v > 1.9){
    warning("Illegal v, should be between 1 and 1.9")
  }                            
  
  powered_exponential <- function(tau){
    pt <- exp(-(tau^v))
    return(pt)
    
  }
  return(powered_exponential)
}


get_matern <- function(sigma2){
  if(v < 1 | v > 5){
    warning("Illegal sigma2, should be between 1 and 1.9")
  }
}

# TODO : LEQ or >< in warnings 

# Plotting 1a) 
xx <- seq(from = 0, to = 10, by = 0.01) 
vs <- seq(from = 1, to = 1.9, by = 0.1)
powered_exponential_functions <- lapply(vs, function(v) get_powered_exponential(v))
corrs <- lapply(powered_exponential_functions, function(x) x(xx))
corrs <- lapply(1:length(corrs), function(x) cbind(corrs[[x]], vs[x]))
corrs <- Reduce(rbind, corrs)
corrs <- as.data.frame(corrs)
colnames(corrs) <- c("corr", "v")
corrs$v <- as.factor(corrs$v)
p1 <- ggplot(corrs, aes(x = rep(xx, length(vs)), y = corr,col = v)) + geom_line() + xlab(TeX("$\\tau$")) + ylab(TeX("$\\rho(\\tau)$"))

