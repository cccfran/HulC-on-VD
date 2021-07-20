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
  
  # HulC
  interval <- cbind(lower = apply(beta_mat, 2, min),
                    upper = apply(beta_mat, 2, max))
  
  interval
}


partialling_out_estimator <- function(data, Y_var, A_var, formula) {
  
  # data
  Y <- data[, Y_var]
  A <- data[, A_var]
  X <- model.matrix(formula, data)
  
  # \hat{E}(A|L)
  E_AL_rf <- regression_forest(X, A)
  E_AL_hat <- predict(E_AL_rf)$predictions
  
  # \hat{E}(Y|L)
  E_YL_rf <- regression_forest(X, Y)
  E_YL_hat <- predict(E_YL_rf)$predictions
  
  # partialling out estimator
  coef_hat <- sum((A - E_AL_hat)*(Y - E_YL_hat))/sum((A - E_AL_hat)^2)
  
  return(list(coefficients = coef_hat))
}

# WIP: sd too big
VD <- function(data, Y_var, A_var, formula,  g = I, g_inv = I) {
  
  # AX formula
  formula_AX <- update(formula, paste0(".~", A_var, "+."))
  
  # data
  Y <- data[, Y_var]
  A <- data[, A_var]
  X <- model.matrix(formula, data)
  AX <- model.matrix(formula_AX, data)
  
  # \hat{E}(A|L)
  E_AL_rf <- regression_forest(X, A)
  E_AL_hat <- predict(E_AL_rf)$predictions
  
  # \hat{E}(Y|A,L)
  E_YAL_rf <- regression_forest(AX, Y)
  E_YAL_hat <- predict(E_YAL_rf)$predictions
  
  # \hat{E}(g(\hat{E}[Y|A,L])|L)
  E_g_E_YAL_hat_A1 <- g(predict(E_YAL_rf, cbind(1, X))$predictions)*E_AL_hat
  E_g_E_YAL_hat_A0 <- g(predict(E_YAL_rf, cbind(0, X))$predictions)*(1-E_AL_hat)
  E_g_E_YAL_hat <- E_g_E_YAL_hat_A0 + E_g_E_YAL_hat_A1
  
  # mu(Y, A, L)
  mu_hat <- g_inv(E_YAL_hat)*(Y-E_YAL_hat) + g(E_YAL_hat) - E_g_E_YAL_hat
  
  # OLS
  new_df <- data.frame(mu_hat = mu_hat,
                       Z = A - E_AL_hat)
  fit <- lm(mu_hat ~ Z - 1, new_df)
  sd <- sqrt(sandwich::sandwich(fit))
  
  return(list(coefficients = fit$coef,
              sd = sd))
}
