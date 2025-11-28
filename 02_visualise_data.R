### Run 01_data wrangling first

# load libraries
library(ggpattern)

# Figure 1 - plot annual ringing counts
ggplot(recap_data_trim) +
  geom_bar(aes(x = year),
           colour = "black") +
  theme_bw(base_size = 14) +
  xlab("Year") +
  ylab("Number of fledglings banded") 

# Figure 2 - plot annual ringing counts, split by seen/not seen as an adult
ggplot(recap_data_trim) +
  geom_bar(aes(x = year, fill = obs_adult),
           colour = "black") +
  scale_fill_manual(values = c("white", "darkgrey")) +
  theme_bw(base_size = 14) +
  xlab("Year") +
  ylab("Number of fledglings banded") +
  labs(fill = "Resighted \nafter \nRecruitment") +
  scale_x_discrete(labels = c('0' = 'No',
                              '1' = 'Yes'))

# Figure 3a - faceted plot of plastic mass
ggplot(recap_data_trim) +
  geom_histogram(aes(x = PlasMassTOT), binwidth = 1,
                 colour = "black") +
  facet_wrap(~ obs_adult,
             labeller = labeller(obs_adult = c('0' = 'Not resighted after banding', 
                                               '1' = 'Resighted after recruitment'))) +
  theme_bw(base_size = 14) +
  xlab("Mass of ingested plastic (grams)") +
  ylab("Observed count")

# Figure 3b - faceted by number of recaptures
ggplot(recap_data_trim) +
  geom_histogram(aes(x = PlasMassTOT), binwidth = 1,
                 colour = "black") +
  facet_wrap(~ n_recap,
             labeller = labeller(n_recap = c('0' = 'Not resighted after banding',
                                             '1' = 'Resighted once after recruitment',
                                             '2' = 'Resighted twice after recruitment',
                                             '3' = 'Resighted 3 times after recruitment',
                                             '4' = 'Resighted 4 times after recruitment',
                                             '5' = 'Resighted 5 times after recruitment'))) +
  theme_bw(base_size = 14) +
  xlab("Mass of ingested plastic (grams)") +
  ylab("Observed count")

# Figure 4 - plot plastic mass by ringing cohorts, facet by seen/not seen as an adult
ggplot(recap_data_trim) +
  geom_jitter(aes(x = year, y = PlasMassTOT), 
              height = 0, width = 0.3) +
  facet_wrap(~ obs_adult,
             labeller = labeller(obs_adult = c('0' = 'Not resighted after banding',
                                               '1' = 'Resighted after recruitment'))) +
  theme_bw(base_size = 14) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  xlab("Year") +
  ylab("Mass of ingested plastic (grams)")




# additional plots

# plot annual ringing counts, split by number of times recaptured as an adult
ggplot(recap_data_trim) +
  geom_bar_pattern(aes(x = year, pattern = n_recap),
                   color = "black",
                   pattern_fill = "black",
                   fill = "white",
                   pattern_angle = 45,
                   pattern_density = 0.2,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  theme_bw(base_size = 14) +
  facet_wrap(~ n_recap) +
  labs(fill = "Resighting \nFrequency")

# plot plastic mass
ggplot(recap_data_trim) +
  geom_histogram(aes(x = PlasMassTOT), binwidth = 1)

# plot non-zero plastic mass
ggplot(recap_data_trim |> filter(PlasMassTOT > 0)) +
  geom_histogram(aes(x = PlasMassTOT), binwidth = 1)

# plot plastic mass, split by seen/not seen as an adult
ggplot(recap_data_trim) +
  geom_histogram(aes(x = PlasMassTOT, colour = obs_adult, fill = obs_adult), binwidth = 1)

# plot plastic mass, split by number of times recaptured as an adult
ggplot(recap_data_trim) +
  geom_histogram(aes(x = PlasMassTOT, colour = n_recap, fill = n_recap), binwidth = 1)

# plot plastic mass by ringing cohorts, facet by seen/not seen as an adult, split by number of times recaptured as an adult
ggplot(recap_data_trim) +
  geom_jitter(aes(x = year, y = PlasMassTOT, colour = n_recap), height = 0, width = 0.3) +
  facet_wrap(~ obs_adult,
             labeller = labeller(obs_adult = c('0' = 'Not resighted',
                                               '1' = 'Resighted'))) +
  theme_grey(base_size = 14) +
  labs(colour = 'Resighting \nFrequency') +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ylab("plastic mass")
