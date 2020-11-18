## Simulation Examples used in the paper

# Example 1

model1 <- function(N, p, family, a, b, c){

  active <- 1:8

  Sigma <- diag(p)

  x <- rmvexp(n = N, p = p, rate = rep(1, p), corr = Sigma); x <- x-1

  true_beta  <- rep(0, length = p)

  B <- rbinom(8, 1, 0.4)

  true_beta[1:8] <- a*((-1)^B)*(4*log(N)/sqrt(N) + abs(rnorm(8, 0, 1)))

  if(family == "gaussian") {

    linpred <-  x%*%cbind(true_beta)

    sigma <- 1

    y <- linpred + sigma*rnorm(N, mean = 0, sd = 1)

    SNR <- snr(x, y, active, family)} else if(family == "binomial") {

      true_beta  <- rep(0, length = p)

      true_beta[1:8] <- b*((-1)^B)*(4*log(N)/sqrt(N) + abs(rnorm(8, 0, 1)))

      linpred <-  x%*%cbind(true_beta)

      prob <- exp(linpred)/(1 + exp(linpred))

      runis <- runif(N, 0, 1)

      y <- ifelse(runis < prob, 1, 0)

      SNR <- snr(x, y, active, family)} else if (family == "poisson") {

        true_beta  <- rep(0, length = p)

        true_beta[1:8] <- c*((-1)^B)*(4*log(N)/sqrt(N) + abs(rnorm(8, 0, 1)))

        linpred <-  exp(x%*%cbind(true_beta))

        y <- rpois(N, lambda = linpred)

        y <- ifelse(y >= quantile(y, .90), quantile(y, .90), y); y <- round(y)

        SNR = snr(x, y, active, family)
      }

  return(list(x = x, y = y, active = active, beta = true_beta,
              SNR = SNR, family = family))
}




# Example 2

model2 <- function(N, p, family, a, b, c){

  active <- 1:8

  Sigma <- diag(p)

  x <- MASS::mvrnorm(n = N, rep(0, p), Sigma)

  true_beta  <- rep(0,length = p)

  B <- rbinom(8, 1, 0.4)

  true_beta[1:8] <- a*((-1)^B)*(4*log(N)/sqrt(N) + abs(rnorm(8, 0, 1)))

  if(family == "gaussian") {

    linpred <-  x%*%cbind(true_beta)

    sigma <- 1

    y <- linpred + sigma*rnorm(N, mean = 0, sd = 1)

    SNR <- snr(x, y, active, family)} else if(family == "binomial") {

      true_beta  <- rep(0, length = p)

      true_beta[1:8] <- b*((-1)^B)*(4*log(N)/sqrt(N) + abs(rnorm(8, 0, 1)))

      linpred  <-   x%*%cbind(true_beta)

      prob <- exp(linpred)/(1 + exp(linpred))

      runis  <-  runif(N,  0, 1)

      y <- ifelse(runis < prob, 1, 0)

      SNR <- snr(x, y, active, family)} else if(family == "poisson"){

        true_beta  <- rep(0, length = p)

        true_beta[1:8] <- c*((-1)^B)*(4*log(N)/sqrt(N) + abs(rnorm(8, 0, 1)))

        linpred <-  exp(x%*%cbind(true_beta))

        y <- rpois(N, lambda = linpred)

        y <- ifelse(y>=quantile(y,.90), quantile(y,.90), y); y <- round(y)

        SNR = snr(x,y,active,family)

      }

  return(list(x = x, y = y, active = active, beta = true_beta,
              SNR = SNR, family = family))
}


# Example 3

model3 <- function(N, p, family, a, b, c){

  active <- 1:5

  Sigma <- diag(p)

  Z  <- MASS::mvrnorm(n=N, rep(0, p), Sigma)

  W  <- MASS::mvrnorm(n=N, rep(0, p), Sigma)

  x  <- matrix(NA, nrow = N, ncol = p)

  for(i in 1:dim(x)[1]){

    for (j in 1:dim(x)[2]){

      if( j <= 5){

        x[i, j] = (Z[i, j] + W[i, j])/sqrt(2)
      } else {

        x[i, j] = (Z[i, j] + sum(Z[i, 1:5]))/2

      }
    }
  }

  true_beta  <- rep(0, length = p)

  for (i in 1:5){

    true_beta[i]<- a*2*i
  }


  if(family == "gaussian"){

    linpred <-  x%*%cbind(true_beta)

    sigma <- 1

    y <- linpred + sigma*rnorm(N, mean = 0, sd = 1)

    SNR <- snr(x, y, active, family)} else if(family == "binomial") {

      true_beta  <-  rep(0, length = p)

      for (i in 1:5){

        true_beta[i] <- b*2*i
      }

      linpred <-  x%*%cbind(true_beta)

      prob <- exp(linpred)/(1 + exp(linpred))

      runis <- runif(N, 0, 1)

      y <- ifelse(runis < prob, 1, 0)

      SNR <- snr(x, y, active, family)} else if(family == "poisson"){

        true_beta   <-  rep(0, length = p)

        for (i in 1:5){

          true_beta[i] <-  c*2*i
        }

        linpred  <-  exp(x%*%cbind(true_beta))

        y  <- rpois(N, lambda = linpred)

        y <- ifelse(y >= quantile(y, .90), quantile(y, .90), y); y <- round(y)

        SNR  <- snr(x, y, active, family)
      }

  return(list(x = x, y = y, active = active, beta = true_beta,
              SNR = SNR, family = family))
}


# Example 4

model4 <- function(N, p, family, a, b, c){

  active <- c(1:7, 501:507)

  Sigma <- matrix(NA, nrow = 500, ncol = 500)

  for(i in 1:dim(Sigma)[1]){
    
    for(j in 1:dim(Sigma)[2]){

      if (i ==j){
        
        Sigma[i, j]  <- 1
        
      }   else {

        Sigma[i,j] <- 0.5^abs(i-j)
      }
    }

  }


  z  <- MASS::mvrnorm(n = N, rep(0, 500), Sigma)

  u1 <- runif(200, -2, 2)

  u <- matrix(NA, nrow = N, ncol = 500)

  u[, 1]<- u1

  for(i in 1:dim(u)[1]){

    for(j in 2:dim(u)[2]){

      u[i, j] <- 0.5*u[i, j-1] + 0.5*runif(1, -2, 2)
    }
  }

  x <- cbind(z, u)

  true_beta  <- rep(0, length = p)

  B <- rbinom(14, 1, 0.4)

  true_beta[c(1:7, 501:507)] <- a*((-1)^B)*(4*log(N)/sqrt(N) + abs(rnorm(14, 0, 1)))

  if(family == "gaussian"){

    linpred <-  x%*%cbind(true_beta)

    sigma <- 1

    y <- linpred + sigma*rnorm(N, mean = 0, sd = 1)

    SNR <- snr(x, y, active, family)} else if(family == "binomial"){

      true_beta  <- rep(0, length = p)

      true_beta[c(1:7, 501:507)] <- b*((-1)^B)*(4*log(N)/sqrt(N) + abs(rnorm(14, 0, 1)))

      linpred <-  x%*%cbind(true_beta)

      prob <- exp(linpred)/(1 + exp(linpred))

      runis <- runif(N,0,1)

      y <- ifelse(runis < prob,1,0)

      SNR <- snr(x, y, active, family)} else if(family == "poisson"){

        true_beta  <- rep(0, length = p)

        true_beta[c(1:7, 501:507)] <- c*((-1)^B)*(4*log(N)/sqrt(N) + abs(rnorm(14, 0, 1)))

        linpred <-  exp(x%*%cbind(true_beta))

        y <- rpois(N,lambda=linpred)

        y <- ifelse(y >= quantile(y, .90), quantile(y, .90), y); y <- round(y)

        SNR <- snr(x, y, active, family)

      }

  return(list(x = x,y = y,active = active,beta = true_beta,
              SNR = SNR,family = family))
}


# Example 5

model5 <- function(N, p, family, a, b, c){

  active <- c(1, 2, 100)

  Sigma <- matrix(NA, nrow = p, ncol = p)

  for(i in 1:dim(Sigma)[1]){
    
    for(j in 1:dim(Sigma)[2]){

      if (i ==j){
        
        Sigma[i, j] <- 1
        
      }   else {

        Sigma[i,j] <- 0.9^abs(i-j)
      }
    }

  }

  x <- MASS::mvrnorm(n = N, rep(0, p), Sigma)

  true_beta  <- rep(0, length = p)

  true_beta[c(1, 2, 100)] <- a*c(-0.5, 1, 0.5)

  if(family == "gaussian"){

    linpred <-  x%*%cbind(true_beta)

    y <- linpred + 0.2*rnorm(N, mean = 0, sd = 1)

    SNR <- snr(x, y, active, family)} else if(family == "binomial"){

      true_beta  <- rep(0, length = p)

      true_beta[c(1, 2, 100)] <- b*c(-0.5, 1, 0.5)

      linpred <-  x%*%cbind(true_beta)

      prob <- exp(linpred)/(1 + exp(linpred))

      runis <- runif(N, 0, 1)

      y <- ifelse(runis < prob, 1, 0)

      SNR <- snr(x, y, active, family)} else if (family == "poisson"){

        true_beta  <- rep(0,length=p)

        true_beta[c(1, 2, 100)] <- c*c(-0.5, 1, 0.5)

        linpred  <-  exp(x%*%cbind(true_beta))

        y <- rpois(N, lambda = linpred)

        y <- ifelse(y >= quantile(y, .90), quantile(y, .90), y); y <- round(y)

        SNR <- snr(x, y, active, family)
      }

  return(list(x = x,y = y, active = active, beta = true_beta,
              SNR = SNR, family = family))
}



### Functions used to produce these examples

rmvexp <- function(n, p, rate, corr){

  ## extract parameters, do sanity checks, deal with univariate case

  if(!is.matrix(corr) || !isSymmetric(corr))
    stop("'corr' must be a symmetric matrix")
  D = ncol(corr)

  Dr = length(rate)
  if(Dr > D)
    warning("'rate' longer than width of 'corr', truncating to fit")
  if(Dr != D)
    rate = rep(rate, length.out=D)

  if(D == 1) rexp(n, rate)

  ## generate standard multivariate normal matrix, convert to CDF

  Z = MASS::mvrnorm(n=n, rep(0, p), corr)
  cdf = pnorm(Z)

  ## convert to exp, return

  sapply(1:D, function(d) qexp(cdf[,d], rate[d]))
}



snr=function(x, y, active, family){

  fit1 <- glm(y~-1, family)

  fit <- glm(y ~ x[, active], family)

  dev <- deviance(fit)

  dev1 <- deviance(fit1)

  SNR <- (dev1-dev)/dev

  return(SNR)
}




## Parameter configurations used in the paper

# Example 1, Gaussian model
model1(N=400, p=1000, family="gaussian", a=1, b=1, c=1)
# Example 2, Gaussian model
model2(N=400, p=1000, family="gaussian", a=1, b=1, c=1)
# Example 3, Gaussian model
model3(N=400, p=1000, family="gaussian", a=1, b=1, c=1)
# Example 4, Gaussian model
model4(N=400, p=1000, family="gaussian", a=1, b=1, c=1)
# Example 5, Gaussian model
model5(N=400, p=1000, family="gaussian", a=1, b=1, c=1)



# Example 1, Binomial model
model1(N=400, p=1000, family="binomial", a=1, b=1, c=1)
# Example 2, Binomial model
model2(N=400, p=1000, family="binomial", a=1, b=1.5, c=1)
# Example 3, Binomial model
model3(N=400, p=1000, family="binomial", a=1, b=1, c=1)
# Example 4, Binomial model
model4(N=400, p=1000, family="binomial", a=1, b=1.5, c=1)
# Example 5, Binomial model
model5(N=400, p=1000, family="binomial", a=1, b=3, c=1)



# Example 1, Poisson model
model1(N=400, p=1000, family="poisson", a=1, b=1, c=.3)
# Example 2, Poisson model
model2(N=400, p=1000, family="poisson", a=1, b=1, c=.3)
# Example 3, Poisson model
model3(N=400, p=1000, family="poisson", a=1, b=1, c=.1)
# Example 4, Poisson model
model4(N=400, p=1000, family="poisson", a=1, b=1, c=.3)
# Example 5, Poisson model
model5(N=400, p=1000, family="poisson", a=1, b=1, c=2)





