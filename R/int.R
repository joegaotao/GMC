
#####################################
## 计算sum-product类型的积分
######################################

### default computing GMC(Y|X)
integration <- function(x, X, Y, h, kernel = "normal"){
  n <- NROW(X)
  kernel <- match.arg(kernel, c("normal", "box", "epanech", "biweight", "triweight"))
  int1 <- 0
  int2 <- 0
  for(i in 1:n){
    point <- (x - X[i])/h
    if(kernel == "normal"){
      int1 <- (Y[i] * dnorm(point)) + int1
      int2 <- dnorm(point) + int2
    }
    if(kernel == "epanech"){
      point <- ifelse(abs(point) > 1, 0, point)
      int1 <- (Y[i] * 3/4 * (1 - point^2)) + int1
      int2 <- 3/4 * (1 - point^2) + int2
    }
    if(kernel == "biweight"){
      point <- ifelse(abs(point) > 1, 0, point)
      int1 <- (Y[i] * 15/16 * (1 - point^2)^2) + int1
      int2 <- 15/16 * (1 - point^2)^2 + int2
    }
    if(kernel == "triweight"){
      point <- ifelse(abs(point) > 1, 0, point)
      int1 <- (Y[i] * 35/32 * (1 - point^2)^3) + int1
      int2 <- 35/32 * (1 - point^2)^3 + int2
    }
  }
  int <- int1^2 / (int2 * (n * h))
  return(int)
}


### default computing sigma
integration1 <- function(x, j, X, Y, h, kernel = "normal"){
  n <- NROW(X)
  kernel <- match.arg(kernel, c("normal", "box", "epanech", "biweight", "triweight"))
  int1 <- 0
  int2 <- 0
  for(i in 1:n){
    point <- (x - X[i])/h
    if(kernel == "normal"){
      int1 <- (Y[i] * dnorm(point)) + int1
      int2 <- dnorm(point) + int2
      int3 <- dnorm((x - X[j])/h)
    }
    if(kernel == "epanech"){
      point <- ifelse(abs(point) > 1, 0, point)
      int1 <- (Y[i] * 3/4 * (1 - point^2)) + int1
      int2 <- 3/4 * (1 - point^2) + int2
      int3 <- 3/4 * (1 - ((x - X[j])/h)^2)
    }
    if(kernel == "biweight"){
      point <- ifelse(abs(point) > 1, 0, point)
      int1 <- (Y[i] * 15/16 * (1 - point^2)^2) + int1
      int2 <- 15/16 * (1 - point^2)^2 + int2
    }
    if(kernel == "triweight"){
      point <- ifelse(abs(point) > 1, 0, point)
      int1 <- (Y[i] * 35/32 * (1 - point^2)^3) + int1
      int2 <- 35/32 * (1 - point^2)^3 + int2
    }
  }
  int <- int3*int1 / (h*int2)
  return(int)
}

### default computing sigma
integration2 <- function(x, j, X, Y, h, kernel = "normal"){
  n <- NROW(X)
  kernel <- match.arg(kernel, c("normal", "box", "epanech", "biweight", "triweight"))
  int1 <- 0
  int2 <- 0
  for(i in 1:n){
    point <- (x - X[i])/h
    if(kernel == "normal"){
      int1 <- (Y[i] * dnorm(point)) + int1
      int2 <- dnorm(point) + int2
      int3 <- dnorm((x - X[j])/h)
    }
    if(kernel == "epanech"){
      point <- ifelse(abs(point) > 1, 0, point)
      int1 <- (Y[i] * 3/4 * (1 - point^2)) + int1
      int2 <- 3/4 * (1 - point^2) + int2
      int3 <- 3/4 * (1 - ((x - X[j])/h)^2)
    }
    if(kernel == "biweight"){
      point <- ifelse(abs(point) > 1, 0, point)
      int1 <- (Y[i] * 15/16 * (1 - point^2)^2) + int1
      int2 <- 15/16 * (1 - point^2)^2 + int2
    }
    if(kernel == "triweight"){
      point <- ifelse(abs(point) > 1, 0, point)
      int1 <- (Y[i] * 35/32 * (1 - point^2)^3) + int1
      int2 <- 35/32 * (1 - point^2)^3 + int2
    }
  }
  int <- int3*int1^2 / (h*int2^2)
  return(int)
}
