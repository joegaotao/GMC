######################################################
#####################################################
# generalized measure of correlation GMC(Y|X)
# @response, Y
# @condition, X
# @bandwidth, bandwidth for kernal
# @kernal, choose kernel method
####################################################
library(MASS)
#####################
#### GMC
##################
GMC <- function(X, Y,bandwidth,kernel = "normal"){
  if(NROW(X) != NROW(Y))
     stop("Two variables should have the same length!")
  n <- NROW(X)
  kernel <- match.arg(kernel, c("normal", "box", "epanech", "biweight", "triweight"))
  
  ## default select the best bandwidth
# bandwidth <- hSelect(X, Y, kernel = kernel)
  ### for computing GMC(Y|X)
   intY2X <- integrate(integration, min(2 * range(X)[1],range(X)[1]), max(2 * range(X)[2], range(X)[2]), X = X, Y = Y, h = bandwidth, kernel = kernel)$value

  meanY <- mean(Y)
  SY <- var(Y) 
  ### for computing GMC(X|Y)
   intX2Y <- integrate(integration, min(2 * range(Y)[1],range(Y)[1]), max(2 * range(Y)[2], range(Y)[2]), X = Y, Y = X, h = bandwidth, kernel = kernel)$value

  meanX <- mean(X)
  SX <- var(X) 
  
  if(kernel == "normal")
    varK <- 1/1
  if(kernel == "epanech")
    varK <- 1/5
  if(kernel == "biweight")
    varK <- 1/7
  if(kernel == "triweight")
    varK <- 1/9
  
  ### com
  gmcY2X <- (intY2X - (meanY + 0)^2)/(SY + bandwidth^2 * varK)
  gmcX2Y <- (intX2Y - (meanX + 0)^2)/(SX + bandwidth^2 * varK)
  
  ##################################################
  ##  calculate the estimation according to 
  ##        Thm 1 and Thm 2
  ##############################################
  ## calculate A
  elementY1 <- -2 * meanY/sqrt(SY + 2 * bandwidth^2 * varK) + 2 * (intY2X - meanY^2) * meanY / (SY + 2 * bandwidth^2 * varK)^{3/2}
  elementY2 <- -(intY2X - meanY^2) / (SY + 2 * bandwidth^2 * varK)^2
  elementX1 <- -2 * meanX/sqrt(SX + 2*bandwidth^2 * varK) + 2 * (intX2Y - meanX^2) * meanX / (SX + 2 * bandwidth^2 * varK)^{3/2}
  elementX2 <- -(intX2Y - meanX^2) / (SX + 2 * bandwidth^2 * varK)^2
  A <- matrix(c(1, elementY1, elementY2, rep(0, 6), 1, elementX1, elementX2), nc = 2)
  
  ## calculate Sigma
  SigmaY1=c()
  SigmaY2=c()
  SigmaX1=c()
  SigmaX2=c()
  for(i in 1:n){
    SigmaY1[i] <- integrate(integration1, min(2 * range(X)[1],range(X)[1]), max(2 * range(X)[2],range(X)[2]), j=i , X = X, Y = Y, h = bandwidth, kernel = kernel)$value
    SigmaY2[i] <- integrate(integration2, min(2 * range(X)[1],range(X)[1]), max(2 * range(X)[2],range(X)[2]), j=i , X = X, Y = Y, h = bandwidth, kernel = kernel)$value
    SigmaX1[i] <- integrate(integration1, min(2 * range(Y)[1],range(Y)[1]), max(2 * range(Y)[2],range(Y)[2]), j=i , X = Y, Y = X, h = bandwidth, kernel = kernel)$value
    SigmaX2[i] <- integrate(integration2, min(2 * range(Y)[1],range(Y)[1]), max(2 * range(Y)[2],range(Y)[2]), j=i , X = Y, Y = X, h = bandwidth, kernel = kernel)$value
  }

  Sigma1 <- (2 * Y * SigmaY1 - SigmaY2)/(SY + 2 * bandwidth^2 * varK)
  Sigma2 <- (2 * X * SigmaX1 - SigmaX2)/(SX + 2 * bandwidth^2 * varK)
  Sigma <- cbind(Sigma1, Y/sqrt(SY + 2 * bandwidth^2 * varK), Y^2,
                 Sigma2, X/sqrt(SX + 2 * bandwidth^2 * varK), X^2)
  Sigma <- cov(Sigma)
  
  ## calculate C0
  C0 <- (intY2X - meanY^2) * (1/(SY + 2 * bandwidth^2 * varK) - 1/(SY+bandwidth^2 * varK)) -
        (intX2Y - meanX^2) * (1/(SX + 2 * bandwidth^2 * varK) - 1/(SX+bandwidth^2 * varK))
  
  ## Thm 1 statistic
    tmp1 <- matrix(c(gmcY2X-(intY2X - meanY^2) * (1/(SY + 2 * bandwidth^2 * varK) - 1/(SY+bandwidth^2 * varK)),
                   gmcX2Y-(intX2Y - meanX^2) * (1/(SX + 2 * bandwidth^2 * varK) - 1/(SX+bandwidth^2 * varK))), nr = 1)
  stat1 <- sqrt(n)*sqrtm(solve(t(A) %*% Sigma %*% A)) %*% t(tmp1)
  stat1 <- as.numeric(stat1)
  
  ## Thm 2 statistic
  one <- matrix(c(1, -1), nr = 1)
  stat2 <- sqrt(n)*(one %*% (t(A) %*% Sigma %*% A) %*% t(one))^(-0.5) * (gmcY2X - gmcX2Y - C0)
  stat2 <- as.numeric(stat2)
  
  return(data.frame(gmcY2X = gmcY2X, gmcX2Y = gmcX2Y, stat1 = stat1, stat2 = stat2))
}

