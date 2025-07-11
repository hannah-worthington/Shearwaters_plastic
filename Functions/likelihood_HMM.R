### Likelihood function
# Objective: To evaluate the negative log-likelihood of the HMM given a set of parameter values and capture histories

# Inputs:
# param - model parameters as given to param_unpack
# model - model structure, default all parameters constant
# histories_i - capture histories for all individuals
# n - number of individuals
# occasions_i - number of occasions per individual
# covariate - individual covariate, default NULL 
# covariate2 - individual covariate, default NULL
# covariate3 - individual covariate, default NULL

# Outputs:
# negloglik - value of the negative log likelihood

likelihood_HMM <- function(param, model = c("constant", "constant", "constant"), histories_i, n, occasions_i, covariate = NULL, covariate2 = NULL, covariate3 = NULL)  {
  
  # storage
  loglik_i <- rep(0, n)
  
  # unpack parameters
  HMM_str <- param_unpack(param, model, n, occasions_i, covariate, covariate2, covariate3, res = F)
  
  # loop over individuals
  for (i in 1:n)  {
    
    # initial contribution
    loglik_partial <- HMM_str$pione[[i]]
      
    # loop over occasions
    for (k in 1:(occasions_i[i] - 1))  {
      loglik_partial <- loglik_partial %*% HMM_str$pmat[[i]][ , , (histories_i[[i]][k] + 1), k]
      loglik_partial <- loglik_partial %*% HMM_str$gamma[[i]][ , , k]
    }
      
    # final occasion
    loglik_partial <- loglik_partial %*% HMM_str$pmat[[i]][ , , (histories_i[[i]][occasions_i[i]] + 1), occasions_i[i]]
      
    # individual contribution to likelihood
    loglik_partial <- loglik_partial %*% matrix(1, nrow = 5, ncol = 1)
      
    # store likelihood contribution for cohort
    if (loglik_partial > 0)  {
      loglik_i[i] <- log(loglik_partial)
    } else if (loglik_partial == 0)  {
      loglik_i[i] <- -10000
    }
  }
    
  # evaluate full likelihood
  negloglik <- -sum(loglik_i)
  
  # return
  return(negloglik)
}