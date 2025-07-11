### Run 01_data wrangling first

# load libraries
library(ggeffects)
library(car)

# Model 1
# Intercept only model, no effect of plastic mass or cohort
bern_mod_intercept <- glm(cbind(obs_adult_num, 1 - obs_adult_num) ~ 1, 
                          family = binomial, 
                          data = recap_data_trim)
summary(bern_mod_intercept)

# plot model and overlay data
ggplot(recap_data_trim) +
  geom_jitter(aes(x = PlasMassTOT, y = obs_adult_num, 
                  colour = obs_adult), width = 0, height = 0.1) +
  geom_smooth(aes(x = PlasMassTOT, y = obs_adult_num),
              method = 'glm', 
              formula = y ~ 1,
              method.args = list(family = binomial)) +
  xlab("plastic mass") +
  ylab("resighting probability") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_grey(base_size = 14) +
  guides(colour = 'none')

# fitted values
predict(bern_mod_intercept, newdata = data.frame('PlasMassTOT' = 0), se.fit = TRUE, type = 'response')



# Model 2
# Ringing cohort intercepts, no effect of plastic mass
bern_mod_intercept_cohort <- glm(cbind(obs_adult_num, 1 - obs_adult_num) ~ year, 
                                 family = binomial, 
                                 data = recap_data_trim)
summary(bern_mod_intercept_cohort)
Anova(bern_mod_intercept_cohort)

# plot model and overlay data
ggplot(recap_data_trim) +
  geom_jitter(aes(x = PlasMassTOT, y = obs_adult_num, 
                  colour = as.factor(obs_adult)), width = 0, height = 0.1) +
  geom_smooth(aes(x = PlasMassTOT, y = obs_adult_num, 
                  fill = year),
              method = 'glm',
              formula = y ~ 1,
              method.args = list(family = binomial)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2))
                       
# faceted
ggplot(recap_data_trim) +
  geom_jitter(aes(x = PlasMassTOT, y = obs_adult_num, 
                  colour = as.factor(obs_adult)), width = 0, height = 0.1) +
  geom_smooth(aes(x = PlasMassTOT, y = obs_adult_num, 
                  fill = year),
              method = 'glm',
              formula = y ~ 1,
              method.args = list(family = binomial)) +
  facet_wrap(~ year) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_grey(base_size = 14) +
  xlab("plastic mass") +
  ylab("resighting probability") +
  guides(colour = 'none', fill = 'none')

# fitted values
predict(bern_mod_intercept_cohort, newdata = data.frame(year = as.factor(2010:2019)), type = 'response', se.fit = TRUE)



# Model 3
# Intercept and plastic mass covariate
bern_mod_plastic <- glm(cbind(obs_adult_num, 1 - obs_adult_num) ~ PlasMassTOT, 
                        family = binomial, 
                        data = recap_data_trim)
summary(bern_mod_plastic)
Anova(bern_mod_plastic)

# plot model and overlay data
ggplot(recap_data_trim) +
  geom_jitter(aes(x = PlasMassTOT, y = obs_adult_num, 
                  colour = as.factor(obs_adult)), width = 0, height = 0.1) +
  geom_smooth(aes(x = PlasMassTOT, y = obs_adult_num),
              method = 'glm',
              formula = y ~ x,
              method.args = list(family = binomial)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_grey(base_size = 14) +
  xlab("plastic mass") +
  ylab("resighting probability") +
  guides(colour = 'none')

# fitted values
predict(bern_mod_plastic, newdata = data.frame(PlasMassTOT = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)), type = 'response', se.fit = TRUE)



# Model 4
# ringing cohort intercepts and plastic mass covariate
bern_mod_plastic_cohort <- glm(cbind(obs_adult_num, 1 - obs_adult_num) ~ as.factor(year) + PlasMassTOT,
                               family = binomial, 
                               data = recap_data_trim)
summary(bern_mod_plastic_cohort)
Anova(bern_mod_plastic_cohort)

# plot model and overlay data 
pred_mod_plastic_cohort <- predict_response(bern_mod_plastic_cohort, 
                                            terms = c('PlasMassTOT [all]', 'year'))
ggplot(recap_data_trim) +
  geom_jitter(aes(x = PlasMassTOT, y = obs_adult_num, 
                  colour = as.factor(obs_adult)), width = 0, height = 0.1) +
  geom_line(data = pred_mod_plastic_cohort,
            aes(x = x, y = predicted,
                group = group)) +
  geom_ribbon(data = pred_mod_plastic_cohort, 
              aes(x = x, ymin = conf.low, ymax = conf.high,
                  fill = group), 
              alpha = 0.1) +
  scale_y_continuous(breaks = seq(0, 1, 0.2))

# faceted
ggplot(recap_data_trim) +
  geom_jitter(aes(x = PlasMassTOT, y = obs_adult_num,
                  colour = as.factor(obs_adult)), width = 0, height = 0.1) +
  scale_colour_manual(values = c("darkgrey", "black")) +
  geom_line(data = pred_mod_plastic_cohort,
            aes(x = x, y = predicted,
                group = group),
            size = 1.25) +
  geom_ribbon(data = pred_mod_plastic_cohort,
              aes(x = x, ymin = conf.low, ymax = conf.high,
                  fill = group),
              alpha = 0.6) +
  scale_fill_manual(values = rep("darkgrey", 10)) +
  facet_wrap(~ group) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 14) +
  xlab("Mass of ingested plastic (grams)") +
  ylab("Estimated mature resighting probability") +
  guides(colour = 'none', fill = 'none')

# fitted values, 2010 cohort
predict(bern_mod_plastic_cohort, newdata = data.frame(PlasMassTOT = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10), year = 2010), type = 'response', se.fit = TRUE)
# fitted values, 2016 cohort
predict(bern_mod_plastic_cohort, newdata = data.frame(PlasMassTOT = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10), year = 2016), type = 'response', se.fit = TRUE)



# Model 5
# ringing cohort intercepts, plastic mass covariate and interaction between plastic mass and ringing cohort
bern_mod_interaction <- glm(cbind(obs_adult_num, 1 - obs_adult_num) ~ as.factor(year)*PlasMassTOT, 
                            family = binomial, 
                            data = recap_data_trim)
summary(bern_mod_interaction)
Anova(bern_mod_interaction)

# plot model and overlay data
ggplot(recap_data_trim) +
  geom_jitter(aes(x = PlasMassTOT, y = obs_adult_num, 
                  colour = as.factor(obs_adult)), width = 0, height = 0.1) +
  geom_smooth(aes(x = PlasMassTOT, y = obs_adult_num, 
                  group = year, fill = year),
              method = 'glm',
              formula = y ~ x,
              method.args = list(family = binomial)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2))

# faceted
ggplot(recap_data_trim) +
  geom_jitter(aes(x = PlasMassTOT, y = obs_adult_num, 
                  colour = as.factor(obs_adult)), width = 0, height = 0.1) +
  geom_smooth(aes(x = PlasMassTOT, y = obs_adult_num, 
                  group = year, fill = year),
              method = 'glm', 
              formula = y ~ x,
              method.args = list(family = binomial)) +
  facet_wrap(~ year) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_grey(base_size = 14) +
  xlab("plastic mass") +
  ylab("resighting probability") +
  guides(colour = 'none', fill = 'none')

# fitted values, 2010 cohort
predict(bern_mod_interaction, newdata = data.frame(PlasMassTOT = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10), year = 2010), type = 'response', se.fit = TRUE)
# fitted values, 2016 cohort
predict(bern_mod_interaction, newdata = data.frame(PlasMassTOT = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10), year = 2016), type = 'response', se.fit = TRUE)



# AICs
AICs <- AIC(bern_mod_intercept, 
            bern_mod_intercept_cohort, 
            bern_mod_plastic, 
            bern_mod_plastic_cohort, 
            bern_mod_interaction)

AICs$AIC - min(AICs$AIC)
