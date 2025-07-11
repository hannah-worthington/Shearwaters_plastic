### Bootstrap cohort histories

# Name: bootstrap_data
# Objective: To bootstrap the capture histories maintaining the number of animals per ringing cohort
# Inputs: histories_i - list of capture histories
#         occasions_i - occasions vector, needs bootstrapped to fit bootstrap model
#         covariate - cohort covariate, default NULL
#         covariate2 - cohort covariate, default NULL
#         covariate3 - cohort covariate, default NULL
# Outputs: boot_histories_i - list of bootstrapped individual histories
#          boot_occasions_i - bootstrapped occasions for individuals
#          boot_covariate - bootstrapped cohort covariate
#          boot_covariate2 - bootstrapped cohort covariate (plastic)
#          boot_covariate3 - bootstrapped cohort covariate (fledge)

bootstrap_data <- function(histories_i, occasions_i, covariate = NULL, covariate2 = NULL, covariate3 = NULL, stratify)  {
  
  # storage
  boot_histories_i <- list()
  boot_occasions_i <- NULL
  boot_covariate <- NULL
  boot_covariate2 <- NULL
  boot_covariate3 <- NULL
  add <- 1
  
  # constants
  cohorts <- unique(stratify)
  
  # for each cohort
  for (c in cohorts)  {
    
    # find individuals in cohort
    orig_cohort_i <- which(stratify == c)
    
    # bootstrap individuals
    boot_sample <- sample(orig_cohort_i, length(orig_cohort_i), replace = T)
    
    for (add_i in boot_sample)  {
      
      boot_histories_i[[add]] <- histories_i[[add_i]]
      boot_occasions_i[add] <- occasions_i[add_i]
      boot_covariate[add] <- covariate[add_i]
      boot_covariate2[add] <- covariate2[add_i]
      boot_covariate3[add] <- covariate3[add_i]
      add <- add + 1
    }
    
  }
  
  # return
  return(list('boot_histories_i' = boot_histories_i,
              'boot_occasions_i' = boot_occasions_i,
              'boot_covariate' = boot_covariate,
              'boot_covariate2' = boot_covariate2,
              'boot_covariate3' = boot_covariate3))
}
