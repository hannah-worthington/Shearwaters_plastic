### bootstrap estimates

# Name: bootstrap_fn
# Objective: To generate bootstrap estimates
# Inputs: nboot - number of bootstraps to perform
#         param - starting values for the optimiser
#         histories_i - capture histories
#         model - model structure for the parameters
#         n - number of individuals
#         occasions_i - number of occasions for each individual
#         covariate - individual covariate, default NULL
#         covariate2 - individual covariate, default = NULL
#         covariate3 - individual covariate, default = NULL
#         stratify - covariate to stratify bootstrap by
#         seed - random seed number
# Outputs: list of bootstrap parameter estimates

bootstrap_fn <- function(nboot, param, histories_i, model, n, occasions_i, covariate = NULL, covariate2 = NULL, covariate3 = NULL, stratify, seed)  {
  
  # storage
  boot_intercept <- rep(0, nboot + 1)
  boot_gradient_plastic <- rep(0, nboot + 1)
  boot_gradient_fledge <- rep(0, nboot + 1)
  boot_p <- matrix(0, nrow = nboot + 1, ncol = 5)
  boot_phi <- rep(0, nboot + 1)
    
  # optimise original data
  opt_orig <- nlm(likelihood_HMM, param, model = model, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = covariate, covariate2 = covariate2, covariate3 = covariate3)
  unpack_orig <- param_unpack(opt_orig$estimate, model, n, occasions_i, covariate, covariate2, covariate3, res = T)
  
  # store the original estimates
  if (model[1] == "constant")  {
    boot_intercept[nboot + 1] <- unpack_orig$r_intercept
    boot_p[(nboot + 1), ] <- unpack_orig$p[[370]][2:6]
    boot_phi[nboot + 1] <- unpack_orig$phi
  } else if (model[1] %in% c("covariate_ss", "covariate_cos", "covariate_sco", "covariate_coco"))  {
    boot_intercept[nboot + 1] <- unpack_orig$r_intercept
    boot_gradient_plastic[nboot + 1] <- unpack_orig$r_gradient
    boot_p[(nboot + 1), ] <- unpack_orig$p[[370]][2:6]
    boot_phi[nboot + 1] <- unpack_orig$phi
  } else if (model[1] %in% c("covariate_sss"))  {
    boot_intercept[nboot + 1] <- unpack_orig$r_intercept
    boot_gradient_plastic[nboot + 1] <- unpack_orig$r_gradient_cov2
    boot_gradient_fledge[nboot + 1] <- unpack_orig$r_gradient_cov3
    boot_p[(nboot + 1), ] <- unpack_orig$p[[370]][2:6]
    boot_phi[nboot + 1] <- unpack_orig$phi
  }
  
  # run bootstrap
  set.seed(seed)
  for (b in 1:nboot)  {
    print(paste('bootstrap', b))
    
    # generate data
    boot_data <- bootstrap_data(histories_i, occasions_i, covariate = covariate, covariate2 = covariate2, covariate3 = covariate3, stratify)
    
    # maximise likelihood
    boot_opt <- nlm(likelihood_HMM, param, model, boot_data$boot_histories_i, n, boot_data$boot_occasions_i, boot_data$boot_covariate, boot_data$boot_covariate2, boot_data$boot_covariate3)
    unpack_boot <- param_unpack(boot_opt$estimate, model, n, boot_data$boot_occasions_i, boot_data$boot_covariate, boot_data$boot_covariate2, boot_data$boot_covariate3, res = T)
    
    # store the bootstrap estimates
    if (model[1] == "constant")  {
      boot_intercept[b] <- unpack_boot$r_intercept
      boot_p[b, ] <- unpack_boot$p[[370]][2:6]
      boot_phi[b] <- unpack_boot$phi
    } else if (model[1] %in% c("covariate_ss", "covariate_cos", "covariate_sco", "covariate_coco"))  {
      boot_intercept[b] <- unpack_boot$r_intercept
      boot_gradient_plastic[b] <- unpack_boot$r_gradient
      boot_p[b, ] <- unpack_boot$p[[370]][2:6]
      boot_phi[b] <- unpack_boot$phi
    } else if (model[1] %in% c("covariate_sss"))  {
      boot_intercept[b] <- unpack_boot$r_intercept
      boot_gradient_plastic[b] <- unpack_boot$r_gradient_cov2
      boot_gradient_fledge[b] <- unpack_boot$r_gradient_cov3
      boot_p[b, ] <- unpack_boot$p[[370]][2:6]
      boot_phi[b] <- unpack_boot$phi
    } 
  }
  
  # return
  return(list('boot_intercept' = boot_intercept,
              'boot_plastic' = boot_gradient_plastic,
              'boot_fledge' = boot_gradient_fledge,
              'boot_p' = boot_p,
              'boot_phi' = boot_phi))
}