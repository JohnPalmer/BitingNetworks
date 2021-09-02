# R visualization of epi curves
library(rhdf5)
library(tidyverse)
library(tidybayes)


D = h5read("data/sims/eb=24000.0_hdn=NouBa_hit=3_human_R0=2.072_mdn=fixed_1.0000000000000001e-7_mls=20_mosquito_R0=1.136_n_h=1000_n_m=8000_n_r=1000_n_s=1000_tp=0.1.jld2", "n_human_infections_reps") 
D = t(D) %>% as_tibble()
names(D) = paste0("rep", 1:ncol(D))
D = D %>% mutate(time = 1:nrow(D)) %>% pivot_longer(cols=-time, names_to="rep", values_to="n_infected" )

ggplot(D %>% filter(time<30), aes(x=time, y=n_infected)) + stat_lineribbon(alpha=1) + scale_fill_brewer(palette = "Reds")+ xlab("time") + ylab("infected humans")
ggsave("plots/epicurve_ribbon_eb=24000.0_hdn=LesCo_hit=3_human_R0=1.12_mdn=fixed_1.0000000000000001e-7_mls=20_mosquito_R0=0.682_n_h=1000_n_m=8000_n_r=2000_n_s=1000_tp=0.035.png", width=4, height=4)

ggplot(D %>% filter(time<30), aes(x=time, y=n_infected, group = rep)) + geom_line(alpha=.05) + xlab("time") + ylab("infected humans")
ggsave("plots/epicurve_eb=24000.0_hdn=LesCo_hit=3_human_R0=1.12_mdn=fixed_1.0000000000000001e-7_mls=20_mosquito_R0=0.682_n_h=1000_n_m=8000_n_r=2000_n_s=1000_tp=0.035.png", width=4, height=4)
