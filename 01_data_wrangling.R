# load libraries
library(tidyverse)

# read in data
encounters_full <- read.csv('./Data/FFSH encounters 20241023.csv')

# trim data to those ringed 2010-2019 inclusive
# remove the very large PlasMassTOT
# remove NA Mass values
# add information on number of recaptures per individual
# add binary indicator for recaptured/not as an adult
recap_data_trim <- encounters_full |> 
  filter(year >= 2010) |>
  filter(year <= 2019) |>
  mutate(year = as.factor(year)) |>
  mutate(group = year) |>
  filter(PlasMassTOT <= 10) |>
  filter(!is.na(Mass)) |>
  mutate(n_recap = X2020 + X2021 + X2022 + X2023 + X2024) |>
  mutate(n_recap = as.factor(n_recap)) |>
  mutate(obs_adult_num = ifelse(n_recap == 0, 0, 1)) |>
  mutate(obs_adult = as.factor(obs_adult_num))

# format into histories starting at ringing
histories_i <- list()
for (i in 1:length(recap_data_trim[,1])) {
  ringing <- recap_data_trim$year[i]
  histories_pull <- filter(recap_data_trim, Band == Band[i]) |>
    select(paste('X', ringing, sep = ''):paste('X', 2024, sep = ''))
  histories_i[[i]] <- as.vector(histories_pull, mode = 'numeric')
}

# store constants
n <- length(histories_i)
occasions_i <- rep(0, n)
for (i in 1:n)  {
  occasions_i[i] <- length(histories_i[[i]])
}

# check data
n
recap_data_trim |>
  count(PlasMassTOT > 0)
recap_data_trim |>
  count(n_recap)
recap_data_trim |>
  count(obs_adult)

# store potential covariates
cov_plastic <- recap_data_trim$PlasMassTOT
cov_cohort <- as.numeric(recap_data_trim$year)
n_cohorts <- length(unique(cov_cohort))
cov_mass <- recap_data_trim$Mass
cov_mass_norm <- (cov_mass - mean(cov_mass))/sd(cov_mass)
