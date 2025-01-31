library(sandwich)
library(grf)
source("HulC/R/auxiliary_functions.R")

HulC <- function(data, alpha = 0.05, FUN, ...)  {
  
  # number of splits
  B <- solve_for_B(alpha, Delta = 0, t = 0)
  
  # create random index
  set.seed(1)
  idx <- 1:nrow(data)
  idx_s <- sample(idx, replace = F)
  idx_list <- split(idx_s, ceiling(idx/length(idx)*B))
  
  # estimates
  beta_list <- lapply(idx_list,
                      function(i) FUN(data = data[i,], ...)$coef )
  beta_mat <- matrix(unlist(beta_list), ncol = length(beta_list[[1]]), byrow = T)
  beta_name <- names(FUN(data = data, ...)$coef)
  
  # HulC
  interval <- data.frame(coef = beta_name, 
                         lower = apply(beta_mat, 2, min),
                         upper = apply(beta_mat, 2, max))
  
  interval
}


partialling_out_estimator <- function(data, Y_var, A_var, formula,
                                      SL.library, ci_level = 0.95) {
  
  # data
  Y <- data[, Y_var]
  A <- data[, A_var]
  X <- as.data.frame(model.matrix(formula, data)[,-1])
  
  set.seed(149)
  # \hat{E}(A|L)
  family <- ifelse(dim(table(A)) == 2, "binomial", "gaussian")
  fit_AL <- SuperLearner(A, X, SL.library=SL.library, family = family)
  ps <- predict(fit_AL)$pred
  
  # \hat{E}(Y|L)
  fit_YL <- SuperLearner(Y, X, SL.library=SL.library)
  mu <- predict(fit_YL)$pred
  
  # partialling out estimator
  fit <- lm(I(Y-mu)~-1+I(A-ps))
  interval <- coefci(fit, vcov=sandwich, level = ci_level)
  
  return(list(coefficients = fit$coef,
              interval = interval))
}


VD <- function(data, Y_var, A_var, formula, SL.library, ci_level = 0.95) {
  
  # AX formula
  formula_AX <- update(formula, paste0(".~", A_var, "+."))
  
  # data
  Y <- data[, Y_var]
  A <- data[, A_var]
  X <- as.data.frame(model.matrix(formula, data)[,-1])
  AX <- as.data.frame(model.matrix(formula_AX, data)[,-1])
  
  # \hat{E}(Y|A,L)
  set.seed(149)
  fit_YAX <- SuperLearner(Y, AX, SL.library=SL.library, family = "binomial")
  
  # mu(Y, A, L)
  QA <- predict(fit_YAX)$pred
  weight<-QA*(1-QA)
  scale<-1/weight
  
  # \hat{E}(A|L)
  family <- ifelse(dim(table(A)) == 2, "binomial", "gaussian")
  fit_AL <- SuperLearner(A, X, SL.library=SL.library, family = family)
  ps <- predict(fit_AL)$pred
  
  if(dim(table(A)) == 2) {
    
    # if A is binary variable
    Q1 <- predict(fit_YAX,
                  newdata = data.frame(cbind(firstep=1,X)), 
                  onlySL = T)$pred
    Q0 <- predict(fit_YAX,
                  newdata = data.frame(cbind(firstep=0,X)), 
                  onlySL = T)$pred
    proj<-ps*qlogis(Q1)+(1-ps)*qlogis(Q0)
    
  } else {
    fit_AL <- SuperLearner(qlogis(QA), X, SL.library=SL.library)
    proj <- predict(fit_AL)$pred    
  }

  mu_hat <- scale*(Y-QA) + qlogis(QA) - proj
  
  # OLS
  fit <- lm(mu_hat~-1+I(A-ps))
  interval <- coefci(fit, vcov=sandwich, level = ci_level)
  
  return(list(coefficients = fit$coef,
              interval = interval))
}


SL.gam.3 <- function(..., deg.gam = 3) SL.gam(..., deg.gam = deg.gam)
SL.gam.4 <- function(..., deg.gam = 4) SL.gam(..., deg.gam = deg.gam)
SL.gam.5 <- function(..., deg.gam = 5) SL.gam(..., deg.gam = deg.gam)
SL.gam.6 <- function(..., deg.gam = 6) SL.gam(..., deg.gam = deg.gam)
SL.gam.7 <- function(..., deg.gam = 7) SL.gam(..., deg.gam = deg.gam)
SL.gam.8 <- function(..., deg.gam = 8) SL.gam(..., deg.gam = deg.gam)
SL.gam.9 <- function(..., deg.gam = 9) SL.gam(..., deg.gam = deg.gam)
SL.gam.10 <- function(..., deg.gam = 10) SL.gam(..., deg.gam = deg.gam)
