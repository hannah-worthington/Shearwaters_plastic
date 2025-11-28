### run 01_data_wrangling first

# source functions
source('./Functions/bootstrap_fn.R')
source('./Functions/bootstrap_data.R')
source('./Functions/likelihood_HMM.R')
source('./Functions/param_unpack.R')

# Model 1
# constant return, constant capture, constant retention
param_start_ccc <- c(0, 0, 0)
model_ccc <- c('constant', 'constant', 'constant')
likelihood_HMM(param_start_ccc, model = model_ccc, histories_i = histories_i, n = n, occasions_i = occasions_i)

opt_ccc <- nlm(likelihood_HMM, param_start_ccc, model = model_ccc, histories_i = histories_i, n = n, occasions_i = occasions_i)
res_ccc <- param_unpack(opt_ccc$estimate, model = model_ccc, n = n, occasions_i = occasions_i, res = T)
AIC_ccc <- 2*opt_ccc$minimum + 2*length(opt_ccc$estimate)
est_ccc <- c(res_ccc$r[1], NA, res_ccc$p[[370]][2:6], res_ccc$phi)

# Model 2
# constant return, temporal capture, constant retention
param_start_ctc <- c(0, rep(0, 5), 0)
model_ctc <- c('constant', 'time', 'constant')
likelihood_HMM(param_start_ctc, model = model_ctc, histories_i = histories_i, n = n, occasions_i = occasions_i)

opt_ctc <- nlm(likelihood_HMM, param_start_ctc, model = model_ctc, histories_i = histories_i, n = n, occasions_i = occasions_i)
res_ctc <- param_unpack(opt_ctc$estimate, model = model_ctc, n = n, occasions_i = occasions_i, res = T)
AIC_ctc <- 2*opt_ctc$minimum + 2*length(opt_ctc$estimate)



# cohort differences, no covariate
# Model 3
# cohort return, constant capture, constant retention
param_start_cocc <- c(rep(0, n_cohorts), 0, 0)
model_cocc <- c('cohort', 'constant', 'constant')
likelihood_HMM(param_start_cocc, model = model_cocc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort)

opt_cocc <- nlm(likelihood_HMM, param_start_cocc, model = model_cocc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort)
res_cocc <- param_unpack(opt_cocc$estimate, model = model_cocc, n = n, occasions_i = occasions_i, covariate = cov_cohort, res = T)
AIC_cocc <- 2*opt_cocc$minimum + 2*length(opt_cocc$estimate)

# Model 4
# cohort return, temporal capture, constant retention
param_start_cotc <- c(rep(0, n_cohorts), rep(0, 5), 0)
model_cotc <- c('cohort', 'time', 'constant')
likelihood_HMM(param_start_cotc, model = model_cotc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort)

opt_cotc <- nlm(likelihood_HMM, param_start_cotc, model = model_cotc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort)
res_cotc <- param_unpack(opt_cotc$estimate, model = model_cotc, n = n, occasions_i = occasions_i, covariate = cov_cohort, res = T)
AIC_cotc <- 2*opt_cotc$minimum + 2*length(opt_cotc$estimate)



# investigating plastic mass as a covariate
# Model 5
# logistic return on plastic (shared intercept, shared slope), constant capture, constant retention
param_start_sscc <- c(0, 0, 0, 0)
model_sscc <- c('covariate_ss', 'constant', 'constant')
likelihood_HMM(param_start_sscc, model = model_sscc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_plastic)

opt_sscc <- nlm(likelihood_HMM, param_start_sscc, model = model_sscc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_plastic)
res_sscc <- param_unpack(opt_sscc$estimate, model = model_sscc, n = n, occasions_i = occasions_i, covariate = cov_plastic, res = T)
AIC_sscc <- 2*opt_sscc$minimum + 2*length(opt_sscc$estimate)

# Model 6
# logistic return on plastic (shared intercept, shared slope), temporal capture, constant retention
param_start_sstc <- c(0, 0, rep(0, 5), 0)
model_sstc <- c('covariate_ss', 'time', 'constant')
likelihood_HMM(param_start_sstc, model = model_sstc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_plastic)

opt_sstc <- nlm(likelihood_HMM, param_start_sstc, model = model_sstc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_plastic)
res_sstc <- param_unpack(opt_sstc$estimate, model = model_sstc, n = n, occasions_i = occasions_i, covariate = cov_plastic, res = T)
AIC_sstc <- 2*opt_sstc$minimum + 2*length(opt_sstc$estimate)

# Model 7
# logistic return on plastic (cohort intercept, shared slope), constant capture, constant retention
param_start_coscc <- c(rep(0, n_cohorts), 0, 0, 0)
model_coscc <- c('covariate_cos', 'constant', 'constant')
likelihood_HMM(param_start_coscc, model = model_coscc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)

opt_coscc <- nlm(likelihood_HMM, param_start_coscc, model = model_coscc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)
res_coscc <- param_unpack(opt_coscc$estimate, model = model_coscc, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic, res = T)
AIC_coscc <- 2*opt_coscc$minimum + 2*length(opt_coscc$estimate)

# Model 8
# logistic return on plastic (cohort intercept, shared slope), temporal capture, constant retention
param_start_costc <- c(rep(0, n_cohorts), 0, rep(0, 5), 0)
model_costc <- c('covariate_cos', 'time', 'constant')
likelihood_HMM(param_start_costc, model = model_costc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)

opt_costc <- nlm(likelihood_HMM, param_start_costc, model = model_costc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)
res_costc <- param_unpack(opt_costc$estimate, model = model_costc, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic, res = T)
AIC_costc <- 2*opt_costc$minimum + 2*length(opt_costc$estimate)

# Model 9
# logistic return on plastic (shared intercept, cohort slope), constant capture, constant retention
param_start_scocc <- c(0, rep(0, n_cohorts), 0, 0)
model_scocc <- c('covariate_sco', 'constant', 'constant')
likelihood_HMM(param_start_scocc, model = model_scocc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)

opt_scocc <- nlm(likelihood_HMM, param_start_scocc, model = model_scocc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)
# restart at new values
opt_scocc <- nlm(likelihood_HMM, opt_scocc$estimate, model = model_scocc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)
res_scocc <- param_unpack(opt_scocc$estimate, model = model_scocc, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic, res = T)
AIC_scocc <- 2*opt_scocc$minimum + 2*length(opt_scocc$estimate)

# Model 10
# logistic return on plastic (shared intercept, cohort slope), temporal capture, constant retention
param_start_scotc <- c(0, rep(0, n_cohorts), rep(0, 5), 0)
model_scotc <- c('covariate_sco', 'time', 'constant')
likelihood_HMM(param_start_scotc, model = model_scotc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)

opt_scotc <- nlm(likelihood_HMM, param_start_scotc, model = model_scotc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)
# restart at new values
opt_scotc <- nlm(likelihood_HMM, opt_scotc$estimate, model = model_scotc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)
res_scotc <- param_unpack(opt_scotc$estimate, model = model_scotc, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic, res = T)
AIC_scotc <- 2*opt_scotc$minimum + 2*length(opt_scotc$estimate)

# Model 11
# logistic return on plastic (cohort intercept, cohort slope), constant capture, constant retention
param_start_cococc <- c(rep(0, n_cohorts), rep(0, n_cohorts), 0, 0)
model_cococc <- c('covariate_coco', 'constant', 'constant')
likelihood_HMM(param_start_cococc, model = model_cococc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)

opt_cococc <- nlm(likelihood_HMM, param_start_cococc, model = model_cococc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)
# restart at new values and repeat as needed
opt_cococc <- nlm(likelihood_HMM, opt_cococc$estimate, model = model_cococc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)
res_cococc <- param_unpack(opt_cococc$estimate, model = model_cococc, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic, res = T)
AIC_cococc <- 2*opt_cococc$minimum + 2*length(opt_cococc$estimate)

# Model 12
# logistic return on plastic (cohort intercept, cohort slope), temporal capture, constant retention
param_start_cocotc <- c(rep(0, n_cohorts), rep(0, n_cohorts), rep(0, 5), 0)
model_cocotc <- c('covariate_coco', 'time', 'constant')
likelihood_HMM(param_start_cocotc, model = model_cocotc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)

opt_cocotc <- nlm(likelihood_HMM, param_start_cocotc, model = model_cocotc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)
# restart at new values
opt_cocotc <- nlm(likelihood_HMM, opt_cocotc$estimate, model = model_cocotc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic)
res_cocotc <- param_unpack(opt_cocotc$estimate, model = model_cocotc, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_plastic, res = T)
AIC_cocotc <- 2*opt_cocotc$minimum + 2*length(opt_cocotc$estimate)



# investigating departure mass as a covariate
# Model 13
# logistic return on normalised mass (shared intercept, shared slope), constant capture, constant retention
param_start_sscc2 <- c(0, 0, 0, 0)
model_sscc2 <- c('covariate_ss', 'constant', 'constant')
likelihood_HMM(param_start_sscc2, model = model_sscc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_mass_norm)

opt_sscc2 <- nlm(likelihood_HMM, param_start_sscc2, model = model_sscc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_mass_norm)
res_sscc2 <- param_unpack(opt_sscc2$estimate, model = model_sscc2, n = n, occasions_i = occasions_i, covariate = cov_mass_norm, res = T)
AIC_sscc2 <- 2*opt_sscc2$minimum + 2*length(opt_sscc2$estimate)

# Model 14
# logistic return on normalised mass (shared intercept, shared slope), temporal capture, constant retention
param_start_sstc2 <- c(0, 0, rep(0, 5), 0)
model_sstc2 <- c('covariate_ss', 'time', 'constant')
likelihood_HMM(param_start_sstc2, model = model_sstc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_mass_norm)

opt_sstc2 <- nlm(likelihood_HMM, param_start_sstc2, model = model_sstc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_mass_norm)
res_sstc2 <- param_unpack(opt_sstc2$estimate, model = model_sstc2, n = n, occasions_i = occasions_i, covariate = cov_mass_norm, res = T)
AIC_sstc2 <- 2*opt_sstc2$minimum + 2*length(opt_sstc2$estimate)

# Model 15
# logistic return on normalised mass (cohort intercept, shared slope), constant capture, constant retention
param_start_coscc2 <- c(rep(0, n_cohorts), 0, 0, 0)
model_coscc2 <- c('covariate_cos', 'constant', 'constant')
likelihood_HMM(param_start_coscc2, model = model_coscc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm)

opt_coscc2 <- nlm(likelihood_HMM, param_start_coscc2, model = model_coscc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm)
res_coscc2 <- param_unpack(opt_coscc2$estimate, model = model_coscc2, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm, res = T)
AIC_coscc2 <- 2*opt_coscc2$minimum + 2*length(opt_coscc2$estimate)

# Model 16
# logistic return on normalised mass (cohort intercept, shared slope), temporal capture, constant retention
param_start_costc2 <- c(rep(0, n_cohorts), 0, rep(0, 5), 0)
model_costc2 <- c('covariate_cos', 'time', 'constant')
likelihood_HMM(param_start_costc2, model = model_costc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm)

opt_costc2 <- nlm(likelihood_HMM, param_start_costc2, model = model_costc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm)
res_costc2 <- param_unpack(opt_costc2$estimate, model = model_costc2, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm, res = T)
AIC_costc2 <- 2*opt_costc2$minimum + 2*length(opt_costc2$estimate)

# Model 17
# logistic return on normalised mass (shared intercept, cohort slope), constant capture, constant retention
param_start_scocc2 <- c(0, rep(0, n_cohorts), 0, 0)
model_scocc2 <- c('covariate_sco', 'constant', 'constant')
likelihood_HMM(param_start_scocc2, model = model_scocc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm)

opt_scocc2 <- nlm(likelihood_HMM, param_start_scocc2, model = model_scocc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm)
res_scocc2 <- param_unpack(opt_scocc2$estimate, model = model_scocc2, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm, res = T)
AIC_scocc2 <- 2*opt_scocc2$minimum + 2*length(opt_scocc2$estimate)

# Model 18
# logistic return on normalised mass (shared intercept, cohort slope), temporal capture, constant retention
param_start_scotc2 <- c(0, rep(0, n_cohorts), rep(0, 5), 0)
model_scotc2 <- c('covariate_sco', 'time', 'constant')
likelihood_HMM(param_start_scotc2, model = model_scotc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm)

opt_scotc2 <- nlm(likelihood_HMM, param_start_scotc2, model = model_scotc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm)
res_scotc2 <- param_unpack(opt_scotc2$estimate, model = model_scotc2, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm, res = T)
AIC_scotc2 <- 2*opt_scotc2$minimum + 2*length(opt_scotc2$estimate)

# Model 19
# logistic return on normalised mass (cohort intercept, cohort slope), constant capture, constant retention
param_start_cococc2 <- c(rep(0, n_cohorts), rep(0, n_cohorts), 0, 0)
model_cococc2 <- c('covariate_coco', 'constant', 'constant')
likelihood_HMM(param_start_cococc2, model = model_cococc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm)

opt_cococc2 <- nlm(likelihood_HMM, param_start_cococc2, model = model_cococc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm)
res_cococc2 <- param_unpack(opt_cococc2$estimate, model = model_cococc2, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm, res = T)
AIC_cococc2 <- 2*opt_cococc2$minimum + 2*length(opt_cococc2$estimate)

# Model 20
# logistic return on normalised mass (cohort intercept, cohort slope), temporal capture, constant retention
param_start_cocotc2 <- c(rep(0, n_cohorts), rep(0, n_cohorts), rep(0, 5), 0)
model_cocotc2 <- c('covariate_coco', 'time', 'constant')
likelihood_HMM(param_start_cocotc2, model = model_cocotc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm)

opt_cocotc2 <- nlm(likelihood_HMM, param_start_cocotc2, model = model_cocotc2, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm)
res_cocotc2 <- param_unpack(opt_cocotc2$estimate, model = model_cocotc2, n = n, occasions_i = occasions_i, covariate = cov_cohort, covariate2 = cov_mass_norm, res = T)
AIC_cocotc2 <- 2*opt_cocotc2$minimum + 2*length(opt_cocotc2$estimate)



# investigate the combination of plastic and mass
# Model 21
# logistic regression on plastic and normalised mass (shared intercept, shared gradients), constant capture, constant retention
param_start_ssscc <- c(0, 0, 0, 0, 0)
model_ssscc <- c('covariate_sss', 'constant', 'constant')
likelihood_HMM(param_start_ssscc, model = model_ssscc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate2 = cov_plastic, covariate3 = cov_mass_norm)

opt_ssscc <- nlm(likelihood_HMM, param_start_ssscc, model = model_ssscc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate2 = cov_plastic, covariate3 = cov_mass_norm)
res_ssscc <- param_unpack(opt_ssscc$estimate, model = model_ssscc, n = n, occasions_i = occasions_i, covariate2 = cov_plastic, covariate3 = cov_mass_norm, res = T)
AIC_ssscc <- 2*opt_ssscc$minimum + 2*length(opt_ssscc$estimate)

# Model 22
# logistic regression on plastic and normalised mass (shared intercept, shared gradients), temporal capture, constant retention
param_start_ssstc <- c(0, 0, 0, rep(0, 5), 0)
model_ssstc <- c('covariate_sss', 'time', 'constant')
likelihood_HMM(param_start_ssstc, model = model_ssstc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate2 = cov_plastic, covariate3 = cov_mass_norm)

opt_ssstc <- nlm(likelihood_HMM, param_start_ssstc, model = model_ssstc, histories_i = histories_i, n = n, occasions_i = occasions_i, covariate2 = cov_plastic, covariate3 = cov_mass_norm)
res_ssstc <- param_unpack(opt_ssstc$estimate, model = model_ssstc, n = n, occasions_i = occasions_i, covariate2 = cov_plastic, covariate3 = cov_mass_norm, res = T)
AIC_ssstc <- 2*opt_ssstc$minimum + 2*length(opt_ssstc$estimate)



# AICs
AICs <- c(AIC_ccc, AIC_ctc, 
          AIC_cocc, AIC_cotc,
          AIC_sscc, AIC_sstc, 
          AIC_coscc, AIC_costc,
          AIC_scocc, AIC_scotc,
          AIC_cococc, AIC_cocotc,
          AIC_sscc2, AIC_sstc2,
          AIC_coscc2, AIC_costc2,
          AIC_scocc2, AIC_scotc2,
          AIC_cococc2, AIC_cocotc2,
          AIC_ssscc, AIC_ssstc)

deltaAICs <- AICs - min(AICs)

# collect results
ModelRes <- cbind('Model' = 1:22,
                  'Intercept' = c(rep(c('shared', 'shared', 'cohort', 'cohort'), times = 5), 'shared', 'shared'),
                  'Gradient' = c(rep(c('-', 'Shared (Plastic Mass)', 'Cohort (Plastic Mass)', 'Shared (Fledging Mass)', 'Cohort (Fledging Mass)'), each = 4), rep('Shared (Plastic Mass) + Shared (Fledging Mass)', times = 2)),
                  'Capture' = rep(c('Constant', 'Temporal'), times = 11),
                  'AIC' = round(AICs, 1),
                  'deltaAIC' = round(deltaAICs, 1))

ModelRes


# model average parameters
res_ctc
res_sstc
res_ssstc

AIC_chosen <- c(AIC_ctc, AIC_sstc, AIC_ssstc)
delta_AIC_chosen <- AIC_chosen - min(AIC_chosen)
relative <- exp(-0.5*delta_AIC_chosen)
weights <- relative/sum(relative)

# weighted parameter estimates
# r intercept
sum(c(res_ctc$r_intercept, res_sstc$r_intercept, res_ssstc$r_intercept)*weights)
# r gradient plastic mass
sum(c(res_sstc$r_gradient, res_ssstc$r_gradient_cov2)*(weights[2:3]/sum(weights[2:3])))
# p
sum(c(res_ctc$p[[370]][2], res_sstc$p[[370]][2], res_ssstc$p[[370]][2])*weights)
sum(c(res_ctc$p[[370]][3], res_sstc$p[[370]][3], res_ssstc$p[[370]][3])*weights)
sum(c(res_ctc$p[[370]][4], res_sstc$p[[370]][4], res_ssstc$p[[370]][4])*weights)
sum(c(res_ctc$p[[370]][5], res_sstc$p[[370]][5], res_ssstc$p[[370]][5])*weights)
sum(c(res_ctc$p[[370]][6], res_sstc$p[[370]][6], res_ssstc$p[[370]][6])*weights)
# phi
sum(c(res_ctc$phi, res_sstc$phi, res_ssstc$phi)*weights)



# bootstrap for standard errors and confidence intervals of three best models
nboot <- 999
rand_seed <- 12345

# Model 2
boot_values_2 <- bootstrap_fn(nboot, opt_ctc$estimate, histories_i, model_ctc, n, occasions_i, stratify = cov_cohort, seed = rand_seed)

SE_r_intercept_2 <- sd(boot_values_2$boot_intercept)
CI_r_intercept_2 <- quantile(boot_values_2$boot_intercept, c(0.025, 0.975))

SE_p_2 <- apply(boot_values_2$boot_p, 2, sd)
CI_p_2 <- apply(boot_values_2$boot_p, 2, quantile, probs = c(0.025, 0.975))

SE_phi_2 <- sd(boot_values_2$boot_phi)
CI_phi_2 <- quantile(boot_values_2$boot_phi, c(0.025, 0.975))

# Model 6
boot_values_6 <- bootstrap_fn(nboot, opt_sstc$estimate, histories_i, model_sstc, n, occasions_i, covariate = cov_plastic, stratify = cov_cohort, seed = rand_seed)

SE_r_intercept_6 <- sd(boot_values_6$boot_intercept)
CI_r_intercept_6 <- quantile(boot_values_6$boot_intercept, c(0.025, 0.975))

SE_r_plastic_6 <- sd(boot_values_6$boot_plastic)
CI_r_plastic_6 <- quantile(boot_values_6$boot_plastic, c(0.025, 0.975))

SE_p_6 <- apply(boot_values_6$boot_p, 2, sd)
CI_p_6 <- apply(boot_values_6$boot_p, 2, quantile, probs = c(0.025, 0.975))

SE_phi_6 <- sd(boot_values_6$boot_phi)
CI_phi_6 <- quantile(boot_values_6$boot_phi, c(0.025, 0.975))

# Model 22
boot_values_22 <- bootstrap_fn(nboot, opt_ssstc$estimate, histories_i, model_ssstc, n, occasions_i, covariate2 = cov_plastic, covariate3 = cov_mass_norm, stratify = cov_cohort, seed = rand_seed)

SE_r_intercept_22 <- sd(boot_values_22$boot_intercept)
CI_r_intercept_22 <- quantile(boot_values_22$boot_intercept, c(0.025, 0.975))

SE_r_plastic_22 <- sd(boot_values_22$boot_plastic)
CI_r_plastic_22 <- quantile(boot_values_22$boot_plastic, c(0.025, 0.975))

SE_r_fledge_22 <- sd(boot_values_22$boot_fledge)
CI_r_fledge_22 <- quantile(boot_values_22$boot_fledge, c(0.025, 0.975))

SE_p_22 <- apply(boot_values_22$boot_p, 2, sd)
CI_p_22 <- apply(boot_values_22$boot_p, 2, quantile, probs = c(0.025, 0.975))

SE_phi_22 <- sd(boot_values_22$boot_phi)
CI_phi_22 <- quantile(boot_values_22$boot_phi, c(0.025, 0.975))


# weighted standard error estimates
# r intercept
sqrt(sum(c(SE_r_intercept_2^2, SE_r_intercept_6^2, SE_r_intercept_22^2)*weights))
# r gradient plastic mass
sqrt(sum(c(SE_r_plastic_6^2, SE_r_plastic_22^2)*(weights[2:3]/sum(weights[2:3]))))
# p
sqrt(sum(c(SE_p_2[1]^2, SE_p_6[1]^2, SE_p_22[1]^2)*weights))
sqrt(sum(c(SE_p_2[2]^2, SE_p_6[2]^2, SE_p_22[2]^2)*weights))
sqrt(sum(c(SE_p_2[3]^2, SE_p_6[3]^2, SE_p_22[3]^2)*weights))
sqrt(sum(c(SE_p_2[4]^2, SE_p_6[4]^2, SE_p_22[4]^2)*weights))
sqrt(sum(c(SE_p_2[5]^2, SE_p_6[5]^2, SE_p_22[5]^2)*weights))
# phi
sqrt(sum(c(SE_phi_2^2, SE_phi_6^2, SE_phi_22^2)*weights))


# weighted confidence intervals
weighted_r_intercept <- rep(0, nboot)
weighted_r_plastic <- rep(0, nboot)
weighted_p <- matrix(0, nboot, 5)
weighted_phi <- rep(0, nboot)
for (b in 1:nboot)  {
  # r intercept
  weighted_r_intercept[b] <- sum(c(boot_values_2$boot_intercept[b], boot_values_6$boot_intercept[b], boot_values_22$boot_intercept[b])*weights)
  # r gradient plastic mass
  weighted_r_plastic[b] <- sum(c(boot_values_6$boot_plastic[b], boot_values_22$boot_plastic[b])*(weights[2:3]/sum(weights[2:3])))
  # p
  weighted_p[b, 1] <- sum(c(boot_values_2$boot_p[b, 1], boot_values_6$boot_p[b, 1], boot_values_22$boot_p[b, 1])*weights)
  weighted_p[b, 2] <- sum(c(boot_values_2$boot_p[b, 2], boot_values_6$boot_p[b, 2], boot_values_22$boot_p[b, 2])*weights)
  weighted_p[b, 3] <- sum(c(boot_values_2$boot_p[b, 3], boot_values_6$boot_p[b, 3], boot_values_22$boot_p[b, 3])*weights)
  weighted_p[b, 4] <- sum(c(boot_values_2$boot_p[b, 4], boot_values_6$boot_p[b, 4], boot_values_22$boot_p[b, 4])*weights)
  weighted_p[b, 5] <- sum(c(boot_values_2$boot_p[b, 5], boot_values_6$boot_p[b, 5], boot_values_22$boot_p[b, 5])*weights)
  # phi
  weighted_phi[b] <- sum(c(boot_values_2$boot_phi[b], boot_values_6$boot_phi[b], boot_values_22$boot_phi[b])*weights)
}

# find quantiles
quantile(weighted_r_intercept, c(0.025, 0.975))
quantile(weighted_r_plastic, c(0.025, 0.975))
apply(weighted_p, 2, quantile, probs = c(0.025, 0.975))
quantile(weighted_phi, c(0.025, 0.975))

# model averages for plotting
pred_plastic <- seq(0, 10, 0.001)
mean(cov_mass)
pred_mass <- 620
pred_mass_norm <- (620 - mean(cov_mass))/sd(cov_mass)
pred_boot_logitr <- matrix(0, nrow = nboot, ncol = length(pred_plastic))
pred_boot_r <- matrix(0, nrow = nboot, ncol = length(pred_plastic))
# loop over bootstraps
for (b in 1:nboot)  {
  # get average parameter estimates
  avg_r_intercept <- sum(c(boot_values_2$boot_intercept[b], boot_values_6$boot_intercept[b], boot_values_22$boot_intercept[b])*weights)
  avg_r_plastic <- sum(c(boot_values_6$boot_plastic[b], boot_values_22$boot_plastic[b])*(weights[2:3]/sum(weights[2:3])))
  avg_r_mass <- boot_values_22$boot_fledge[b]
  # loop over plastic values
  for (p in 1:length(pred_plastic))  {
    # estimate r
    pred_boot_logitr[b, p] <- avg_r_intercept + avg_r_plastic*pred_plastic[p] + avg_r_mass*pred_mass_norm
    pred_boot_r[b, p] <- 1 / (1 + exp(-pred_boot_logitr[b, p]))
  }
}

# get se for each plastic value
pred_se_logitr <- apply(pred_boot_logitr, 2, sd)
pred_ci_logitr <- apply(pred_boot_logitr, 2, quantile, probs = c(0.025, 0.975))
pred_se_r <- apply(pred_boot_r, 2, sd)
pred_ci_r <- apply(pred_boot_r, 2, quantile, probs = c(0.025, 0.975))
  
# MLE
MLE_avg_intercept <- sum(c(res_ctc$r_intercept, res_sstc$r_intercept, res_ssstc$r_intercept)*weights)
MLE_avg_plastic <- sum(c(res_sstc$r_gradient, res_ssstc$r_gradient_cov2)*(weights[2:3]/sum(weights[2:3])))
MLE_avg_mass <- res_ssstc$r_gradient_cov3

MLE_logitr <- MLE_avg_intercept + MLE_avg_plastic*pred_plastic + MLE_avg_mass*pred_mass_norm
MLE_r <- 1 / (1 + exp(-MLE_logitr))
  
# lower and upper CI
# pred_lower_logitr <- MLE_logitr + qnorm(0.025)*pred_se_logitr
# pred_upper_logitr <- MLE_logitr + qnorm(0.975)*pred_se_logitr
# pred_lower_r <- 1 / (1 + exp(-pred_lower_logitr))
# pred_upper_r <- 1 / (1 + exp(-pred_upper_logitr))

pred_lower_quant_logitr <- apply(pred_boot_logitr, 2, quantile, probs = 0.025)
pred_upper_quant_logitr <- apply(pred_boot_logitr, 2, quantile, probs = 0.975)
pred_lower_quant_r <- 1 / (1 + exp(-pred_lower_quant_logitr))
pred_upper_quant_r <- 1 / (1 + exp(-pred_upper_quant_logitr))

# visualise the estimated probability of return for a range of plastic values assuming an average fledging mass
ggplot() +
  geom_jitter(data = recap_data_trim,
              aes(x = PlasMassTOT, y = obs_adult_num, 
                  colour = as.factor(obs_adult)), width = 0, height = 0.1) +
  scale_colour_manual(values = c("darkgrey", "black")) +
  geom_line(aes(x = pred_plastic, y = MLE_r), colour = 'black', linewidth = 1.25) +
  geom_ribbon(aes(x = pred_plastic, ymin = pred_lower_quant_r, ymax = pred_upper_quant_r), fill = 'darkgrey', alpha = 0.6) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 14) +
  xlab("Mass of ingested plastic (grams)") +
  ylab("Estimated probability of recruitment (r)") +
  guides(colour = 'none')

examples <- c(1, 2, 3, 6, 11, 21, 51, 101, 201, 501, 1001, 2001, 5001, 10001)
pred_plastic[examples]
MLE_r[examples]
pred_se_r[examples]
pred_lower_quant_r[examples]
pred_upper_quant_r[examples]

which(pred_lower_quant_r < 0.01)[1]
pred_plastic[2019]
