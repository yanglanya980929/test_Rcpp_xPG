loglike_2 <- function(xt, params, y){
  vy <- params[1]
  dnorm(y, mean = xt, sd = sqrt(vy), log=TRUE)
}

simPara.a <- function(b, vx, path){
  N <- length(path)
  m <- (sum(c((path - b)[-1],0)*path))/(sum(path**2)-path[N]**2)
  sigma2 <- vx/(sum(path**2)-path[N]**2)
  return(rnorm(1, mean = m, sd = sqrt(sigma2)))
}