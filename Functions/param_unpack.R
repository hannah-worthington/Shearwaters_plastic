### Unpack function
# Objective: To unpack and transform a vector of parameter values for the HMM likelihood

# Inputs:
# param
#   - logit(r), probability of return (logit scale; constant, cohort, covariate_ss, covariate_cos, covariate_sco, covariate_coco, covariate_sss)
#   - logit(p), probability of capture (logit scale; constant, time)
#   - logit(phi), probability of retention (logit scale; constant)
# model - model structure, default all parameters constant
# n - number of individuals
# occasions_i - number of occasions for each individual
# covariate - individual covariate, default = NULL
# covariate2 - individual covariate, default = NULL
# covariate3 - individual covariate, default = NULL
# res
#   - TRUE returns only parameters of interest
#   - FALSE returns parameters and HMM structures

# Outputs:
# res = TRUE or FALSE
#   - r, probability of return, parameters related to r
#   - beta, probability of arrival on each occasion
#   - p, probability of capture
#   - phi, probability of retention
# res = FALSE, above plus
#   - betastar, conditional probabilities of arrival for each occasion
#   - pione, initial state probabilities
#   - gamma, transition probability matrices
#   - pmat, observation process matrices

param_unpack <- function(param, model = c("constant", "constant", "constant"), n, occasions_i, covariate = NULL, covariate2 = NULL, covariate3 = NULL, res = F)  {

  # probability of return
  if (model[1] == "constant")  {  # intercept only model
    logitr <- param[1]
    param <- param[-1]
    r <- rep(1 / (1 + exp(-logitr)), n)
  } else if (model[1] == 'cohort')  {  # ringing cohort
    n_cohorts <- length(unique(covariate))
    logitr_cohort <- param[1:n_cohorts]
    param <- param[-c(1:n_cohorts)]
    r_cohort <- 1 / (1 + exp(-logitr_cohort))
    r <- rep(0, n)
    for (i in 1:n)  {
      r[i] <- r_cohort[covariate[i]] 
    }
  } else if (model[1] == "covariate_ss")  {  # logistic regression shared intercept, shared covariate slope
    r_intercept <- param[1]
    r_gradient <- param[2]
    param <- param[-c(1, 2)]
    r <- 1 / (1 + exp(-(r_intercept + r_gradient*covariate)))
  } else if (model[1] == "covariate_cos")  {  # logistic regression cohort intercept, shared covariate slope
    n_cohorts <- length(unique(covariate))
    r_intercept <- param[1:n_cohorts]
    r_gradient <- param[n_cohorts + 1]
    param <- param[-c(1:(n_cohorts + 1))]
    r <- rep(0, n)
    for (i in 1:n)  {
      r[i] <- 1 / (1 + exp(-(r_intercept[covariate[i]] + r_gradient*covariate2[i])))
    }
  } else if (model[1] == "covariate_sco")  {  # logistic regression shared intercept, cohort covariate slope
    n_cohorts <- length(unique(covariate))
    r_intercept <- param[1]
    r_gradient <- param[2:(n_cohorts + 1)]
    param <- param[-c(1:(n_cohorts + 1))]
    r <- rep(0, n)
    for (i in 1:n)  {
      r[i] <- 1 / (1 + exp(-(r_intercept + r_gradient[covariate[i]]*covariate2[i])))
    }
  } else if (model[1] == 'covariate_coco')  {  # logistic regression cohort intercept and cohort covariate slope
    n_cohorts <- length(unique(covariate))
    r_intercept <- param[1:n_cohorts]
    param <- param[-c(1:n_cohorts)]
    r_gradient <- param[1:n_cohorts]
    param <- param[-c(1:n_cohorts)]
    r <- rep(0, n)
    for (i in 1:n)  {
      r[i] <- 1 / (1 + exp(-(r_intercept[covariate[i]] + r_gradient[covariate[i]]*covariate2[i])))
    }
  } else if (model[1] == 'covariate_sss')  {  # logistic regression with two covariates
    r_intercept <- param[1]
    r_gradient <- param[2]
    r_gradient2 <- param[3]
    param <- param[-c(1:3)]
    r <- rep(0, n)
    for (i in 1:n)  {
      r[i] <- 1 / (1 + exp(-(r_intercept + r_gradient*covariate2[i] + r_gradient2*covariate3[i])))
    }
  }
  
  # arrival probabilities
  beta <- list()
  for (i in 1:n)  {
    beta_i <- rep(0, (occasions_i[i] - 1))
    beta_i[4] <- 1  # arrive with certainty on 4th season after ringing
    beta[[i]] <- beta_i
  }
  
  # conditional arrival probabilities
  betastar <- list()
  for (i in 1:n)  {
    n_beta <- length(beta[[i]])
    betastar[[i]] <- rep(0, n_beta) 
    for (k in 1:n_beta)  {
      if (sum(beta[[i]][k:n_beta]) > 0)  {
        betastar[[i]][k] <- beta[[i]][k] / sum(beta[[i]][k:n_beta])
      }
    }
  }
  
  # capture probabilities
  if (model[2] == "constant")  {
    logitp <- param[1]
    param <- param[-1]
    pnonzero <- 1 / (1 + exp(-logitp))
    p <- list()
    for (i in 1:n)  {
      p_i <- rep(0, occasions_i[i])
      p_i[1] <- 1
      p_i[(occasions_i[i] - 4):occasions_i[i]] <- pnonzero # specific to number of recapture occasions
      p[[i]] <- p_i
    }
  } else if (model[2] == "time")  {
    logitp <- param[1:5]
    param <- param[-c(1:5)]
    pnonzero_time <- 1 / (1 + exp(-logitp))
    p <- list()
    for (i in 1:n)  {
      p_i <- rep(0, occasions_i[i])
      p_i[1] <- 1
      p_i[(occasions_i[i] - 4):occasions_i[i]] <- pnonzero_time 
      p[[i]] <- p_i
    }
  } 
  
  # retention probabilities
  logitphi <- param[1]
  phi <- 1 / (1 + exp(-logitphi)) # constant
  
  # state probabilities at time 0
  pione <- list()
  for (i in 1:n)  {
    pione[[i]] <- rep(0, 5)
    pione[[i]][1] <- 1
  }
  
  # transition probability matrices
  gamma <- list()
  for (i in 1:n)  {
    gamma[[i]] <- array(0, dim = c(5, 5, (occasions_i[i] - 1)))
    gamma[[i]][1, 2, ] <- r[i]
    gamma[[i]][1, 5, ] <- 1 - r[i]
    gamma[[i]][2, 2, ] <- 1 - betastar[[i]]
    gamma[[i]][2, 3, ] <- betastar[[i]]
    gamma[[i]][3, 3, ] <- phi
    gamma[[i]][3, 4, ] <- 1 - phi
    gamma[[i]][4, 4, ] <- 1
    gamma[[i]][5, 5, ] <- 1
  }
  
  # observation matrices
  pmat <- list()
  for (i in 1:n)  {
    pmat[[i]] <- array(0, dim = c(5, 5, 2, occasions_i[i]))
    pmat[[i]][2, 2, 1, ] <- 1
    pmat[[i]][3, 3, 1, ] <- 1 - p[[i]]
    pmat[[i]][4, 4, 1, ] <- 1
    pmat[[i]][5, 5, 1, ] <- 1
    pmat[[i]][1, 1, 2, ] <- 1
    pmat[[i]][3, 3, 2, ] <- p[[i]]
  }
  
  # return all parameters and structures if res is F
  if (res == F)  {
    return(list('r' = r,
                'beta' = beta,
                'betastar' = betastar, 
                'p' = p, 
                'phi' = phi, 
                'pione' = pione,
                'gamma' = gamma,
                'pmat' = pmat))
  
  # only return parameter values if res is T
  } else if (res == T)  {
    if (model[1] == "constant")  {
      return(list('r_intercept' = logitr,
                  'r' = r,
                  'beta' = beta,
                  'p' = p, 
                  'phi' = phi))
    } else if (model[1] == "cohort")  {
      return(list('r_intercept' = logitr_cohort,
                  'r' = r_cohort,
                  'beta' = beta,
                  'p' = p,
                  'phi' = phi))
    } else if (model[1] %in% c("covariate_ss", "covariate_cos", "covariate_sco", "covariate_coco"))  {
      return(list('r_intercept' = r_intercept,
                  'r_gradient' = r_gradient,
                  'r' = r,
                  'beta' = beta,
                  'p' = p,
                  'phi' = phi))
    } else if (model[1] %in% c("covariate_sss"))  {
      return(list('r_intercept' = r_intercept,
                  'r_gradient_cov2' = r_gradient,
                  'r_gradient_cov3' = r_gradient2,
                  'r' = r,
                  'beta' = beta,
                  'p' = p,
                  'phi' = phi))
    }
  }
}