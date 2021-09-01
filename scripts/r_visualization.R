# Visualization of sim results in R

# dependencies ####
library(tidyverse)
library(RColorBrewer)
library(tidybayes)


plot_data <- read_csv("data/sim_summaries/combo_plot_AR_R0_expected_bites=4000.0_human_infection_time=3_mosquito_life_span=20_n_humans=1000_n_mosquitoes=4000_n_reps=1000_n_steps=1000_this_set_name=bcn_probs_transmission_prob=0.15.csv")

ggplot(plot_data, aes(x=value, y=variable)) + stat_halfeye() + facet_wrap(~measure, scale="free_x")

pop_dist = read_csv("../BarcelonaTiger/data/proc/biting_networks_veri_census_sections_bcn_pop_props.csv") %>% rename(variable = District, value=vri) %>% mutate(measure="VRI") %>% select(variable, value, measure)

pop_dist_all = pop_dist %>% mutate(variable = "Barcelona (all)")

pop_dist = bind_rows(pop_dist, pop_dist_all)

ggplot(bind_rows(plot_data, pop_dist), aes(x=value, y=variable)) + stat_halfeye() + facet_wrap(~measure, scale="free_x")


hbd = read_csv("data/sim_summaries/human_bite_distributions_expected_bites=4000.0_human_infection_time=3_mosquito_life_span=20_n_humans=1000_n_mosquitoes=4000_n_reps=1000_n_steps=1000_this_set_name=bcn_probs_transmission_prob=0.15.csv") %>% pivot_longer(cols=everything(), names_to="variable", values_to="value") %>% mutate(measure = "Bites per person")

ggplot(hbd, aes(x=value, y=variable)) + stat_halfeye() #+ facet_wrap(~measure, scale="free_x"))

plot_data = bind_rows(hbd, plot_data, pop_dist) %>% mutate(measure_f = factor(measure, levels=c("VRI", "Bites per person", "R0", "AR")))

these_medians = plot_data %>% group_by(measure_f) %>% summarise(med = median(value))

ggplot(plot_data, aes(x=value, y=variable, fill=measure_f)) + stat_halfeye(normalize="groups", alpha=.7, point_interval = mean_qi,.width = c(0.66, 0.95)) + facet_wrap(~measure_f, scale="free_x", nrow = 1) + theme(legend.position="none") + scale_fill_brewer(palette = "Set1") + geom_vline(data = these_medians, aes(xintercept=med), lwd=.5, color="yellow") 
