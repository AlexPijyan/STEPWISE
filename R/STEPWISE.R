#' A variable selection and parameter estimation method
#'
#' This function implements STEPWISE method that performs a variable selection and consistent estimation of
#' the model parameters in Generalized Linear Models (GLM).
#' @param x A set of predictor variables. Should be either a numeric matrix or a data frame.
#' @param y The outcome.
#' @param family A description of the error distribution and link function to be used in the model.
#' This can be a character string naming a family function, a family function or the result of a call to a family function.
#' @param eta1 A stopping criteria for controlling the false negatives and positives at the forward step of the method.
#'  Can be any non-negative number. By default, eta1 = 0.
#' @param eta2 A stopping criteria for controlling the false negatives and positives at the backward step of the method.
#'  Can be any non-negative number. By default, eta2 = 1.
#' @param initvars An initial set of predictors included in the model. Should be a character vector with variable names.
#'  By default, initvars = NULL.
#' @return A list of four elements: "finalset", "finalset2", "coef1" and "coef2". "Finalset" is a set of predictors selected
#'  by the forward step of the STEPWISE method; "Finalset2" is a set of predictors selected after completing both forward and
#'  backward steps of the STEPWISE method and should be considered as the final result of the method. "Coef1" and "Coef2"
#'  are parameter estimates of the GLM models built with "Finalset" and "Finalset2", respectively.
#' @example
#'  STEPWISE(x, y, family = "gaussian", eta1 = 1, eta2 = 3.5, initvars = c("Age", "Gender"))
#'



STEPWISE <- function(x, y, family, eta1 = 0, eta2 = 1, initvars = NULL){
  
  nc <- parallel::detectCores(logical = F)
  
  cl <- parallel::makeCluster(nc-1, type = "SOCK")
  
  doParallel::registerDoParallel(cl)
  
  `%dopar%` <- foreach::`%dopar%`
  
  `%>%` <- magrittr::`%>%`
  
  N <- x %>% nrow(); p <- x %>% ncol(); eta1 <- eta1; eta2 <- eta2
  
  EBIC <- min(N,p) %>% numeric(); S <- min(N,p) %>% numeric(); x <- x %>% as.matrix()
  
  
  
  if(initvars %>% is.null()) {
    
    M <- foreach::foreach (i = c(1:dim(x)[2]), .combine = c, .packages = c("MASS", "magrittr")) %dopar% {
      
      glm(y ~ x[, i], family = family) %>% logLik()
      
    }
    
    S[1] <- which.max(M); newset <- S[1]
    
    fit1 <- glm(y ~ x[, S[1]], family = family)
    
    MM <- fit1 %>% logLik()
    
    EBIC[1] <- -2*(MM) + log(N) + eta1*log(p)
    
    
    for(k in 2:min(N,p)) {
      
      M <- foreach::foreach (i = c(1:dim(x)[2]), .combine = c, .packages = c("MASS", "magrittr")) %dopar% {
        
        glm(y ~ x[,i] + as.matrix(x[, newset]), family=family) %>% logLik()
        
      }
      
      M[newset] <- -10^10; S[k] <- which.max(M); newset <- S[1:k]
      
      fit1 <- glm(y ~ as.matrix(x[, newset]), family = family)
      
      MM <- fit1 %>% logLik()
      
      EBIC[k] <- -2*(MM) + k*log(N) + k*eta1*log(p)
      
      if(coef(fit1)[-1] %>% abs() %>% sum() > 30)  break 
      
      if(fit1$converged == FALSE | EBIC[k]>EBIC[k-1] | k == min(N,p))  break 
      
    }
    
    
    finalset <- S[1:(k-1)]
    
    x_new <- x[, finalset, drop = FALSE]
    
    p2 <- x_new %>% ncol(); EBIC2 <- p2 %>% numeric()
    
    fit2 <- glm(y ~ x_new, family = family)
    
    MM <- fit2 %>% logLik()
    
    EBIC2[1] <- -2*(MM) + p2*eta2*log(N)
    
    
    for (j in 2:p2){
      
      if(p2 == 1) break 
      
      M <- dim(x_new)[2] %>% numeric()
      
      for( i in 1:dim(x_new)[2]){
        
        M[i] <- glm(y ~ as.matrix(x_new[, -i, drop=FALSE]), family = family) %>% logLik()
        
      }
      
      
      fit2 <- glm(y ~ as.matrix(x_new[, - which.max(M)]), family = family)
      
      MM <- fit2 %>% logLik()
      
      EBIC2[j] <- -2*(MM) + (p2-j+1)*eta2*log(N)
      
      if(fit2$converged == FALSE | EBIC2[j] > EBIC2[j-1] | j == p2)  break 
      
      x_new <- x_new[, - which.max(M), drop = FALSE]
      
    }
    
    
    finalset2 <- x_new %>% colnames()
    
    fit1 <- glm(y ~ as.matrix(x[, finalset]), family = family) %>% summary()
    
    fit2 <- glm(y ~ as.matrix(x[, finalset2]), family = family) %>% summary()
    
    
  } else {
    
    M <- foreach::foreach (i=c(1:dim(x)[2]), .combine = c, .packages = c("MASS", "magrittr")) %dopar% {
      
      glm(y ~ as.matrix(x[, i])+ as.matrix(x[, initvars]), family = family) %>% logLik()
      
    }
    
    S[1] <- which.max(M);  newset <- S[1]
    
    fit1 <- glm(y ~ as.matrix(x[, S[1]])+ as.matrix(x[, initvars]), family = family)
    
    MM <- fit1 %>% logLik()
    
    EBIC[1] = -2*(MM) + log(N) + eta1*log(p)
    
    
    for(k in 2:min(N,p)){
      
      M <- foreach::foreach (i = c(1:dim(x)[2]), .combine = c, .packages = c("MASS", "magrittr")) %dopar% {
        
        glm(y ~ as.matrix(x[,i]) + as.matrix(x[, newset])+ as.matrix(x[, initvars]), family = family) %>% logLik()
        
      }
      
      M[newset] <- -10^10; S[k] <- which.max(M); newset <- S[1:k]
      
      fit1 <- glm(y ~as.matrix(x[, newset])+ as.matrix(x[, initvars]), family = family)
      
      MM <- fit1 %>% logLik()
      
      EBIC[k] <- -2*(MM) + k*log(N) + k*eta1*log(p)
      
      
      if(coef(fit1)[-1] %>% abs() %>% sum() > 30)  break 
      
      if(fit1$converged == FALSE | EBIC[k]>EBIC[k-1] | k == min(N,p)) break
      
    }
    
    finalset <- S[1:(k-1)]
    
    x_new <- x[, finalset, drop = FALSE]
    
    p2 <- x_new %>% ncol(); EBIC2 <- p2 %>% numeric()
    
    fit2 <- glm(y ~ as.matrix(x_new) + as.matrix(x[, initvars]), family = family)
    
    MM <- fit2 %>% logLik()
    
    EBIC2[1] <- -2*(MM) + p2*eta2*log(N)
    
    
    for (j in 2:p2){
      
      if(p2 == 1) break 
      
      M <- dim(x_new)[2] %>% numeric()
      
      for( i in 1:dim(x_new)[2]){
        
        M[i] <- glm(y ~ as.matrix(x_new[, -i, drop=FALSE]) + as.matrix(x[, initvars]), family = family) %>% logLik()
        
      }
      
      fit2 <- glm(y ~ as.matrix(x_new[, - which.max(M)])+ as.matrix(x[, initvars]), family = family)
      
      MM <- fit2 %>% logLik()
      
      EBIC2[j] <- -2*(MM)+(p2-j+1)*eta2*log(N)
      
      if(fit2$converged == FALSE | EBIC2[j] > EBIC2[j-1] | j == p2) break 
      
      x_new <- x_new[, - which.max(M), drop = FALSE]
      
    }
    
    
    finalset2 <- x_new %>% colnames()
    
    fit1 <- glm(y ~ as.matrix(x[, finalset]), family = family) %>% summary()
    
    fit2 <- glm(y ~ as.matrix(x[, finalset2]), family = family) %>% summary()
    
    
  }
  
  
  
  list(finalset = finalset, finalset2 = finalset2, coef1 = fit1 %>% coefficients(), coef2 = fit2 %>% coefficients())
  
}

