# Visualization of sim results in R

# dependencies ####
library(tidyverse)
library(RColorBrewer)

plot_data <- read_csv("data/sim_summaries/combo_plot_AR_R0_expected_bites=4000.0_human_infection_time=3_mosquito_life_span=20_n_humans=1000_n_mosquitoes=4000_n_reps=1000_n_steps=1000_this_set_name=bcn_probs_transmission_prob=0.15csv")

