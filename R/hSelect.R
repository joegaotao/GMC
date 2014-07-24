
####################
### choice of bandwidth h
## using ks package
###################
hSelect <- function(X, Y, kernel = "normal"){
  n <- NROW(X)
  kernel <- match.arg(kernel, c("normal", "box", "epanech", "biweight", "triweight"))
  value <- optimize(FH, c(0, 10), X, Y, kernel=kernel)
  h.optim <- value$minimum
  return(h.optim)
}

FH <- function(h, X, Y, kernel){
  kernel <- match.arg(kernel, c("normal", "box", "epanech", "biweight", "triweight"))
  w1 <- sqrt(1/(var(Y)))
  w2 <- sqrt(1/(var(X)))
  sumY <- vector()
  sumX <- vector()
  for(i in 1:n){
    point1 <- (X[i] - X[-i])/h
    point2 <- (Y[i] - Y[-i])/h
    if(kernel == "normal"){
      intY <- sum(Y[-i] * dnorm(point1)) / sum(dnorm(point1))
      intX <- sum(X[-i] * dnorm(point2)) / sum(dnorm(point2))
    }
    ### except normal, other kernels are 0 outsie the interval [-1, 1]
    if(kernel == "epanech"){
      point1 <- ifelse(abs(point1) > 1, 0, point1)
      point2 <- ifelse(abs(point2) > 1, 0, point2)
      intY <- sum(Y[-i] * 3/4 * (1 - point1^2)) / sum(3/4 * (1 - point1^2))
      intX <- sum(X[-i] * 3/4 * (1 - point2^2)) / sum(3/4 * (1 - point2^2))
    }
    if(kernel == "biweight"){
      point1 <- ifelse(abs(point1) > 1, 0, point1)
      point2 <- ifelse(abs(point2) > 1, 0, point2)
      intY <- sum(Y[-i] * 15/16 * (1 - point1^2)^2) / sum(15/16 * (1 - point1^2)^2) 
      intX <- sum(X[-i] * 15/16 * (1 - point2^2)^2) / sum(15/16 * (1 - point2^2)^2)
    }
    if(kernel == "triweight"){
      point1 <- ifelse(abs(point1) > 1, 0, point1)
      point2 <- ifelse(abs(point2) > 1, 0, point2)
      intY <- sum(Y[-i] * 35/32 * (1 - point1^2)^3) / sum(35/32 * (1 - point1^2)^3)
      intX <- sum(X[-i] * 35/32 * (1 - point2^2)^3) / sum(35/32 * (1 - point2^2)^3)
    }
    sumY[i] <- (Y[i] - intY)^2
    sumX[i] <- (X[i] - intX)^2
  }
  w1/n * sum(sumY) + w2/n * sum(sumX)
}



