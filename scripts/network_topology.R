# Visualization of sim results in R

rm(list=ls())

# dependencies ####
library(tidyverse)
library(RColorBrewer)
library(tidybayes)

# emporical distributions ####

these_params = "_eb=1600.0_hit=3_mls=20_n_h=100_n_m=800_n_r=10000_n_s=1000_this_set_name=bcn_probs_updated_2022_06_5_tp=0.1"


plot_data <- read_csv(paste0("data/sim_summaries/combo_plot_AR_R0", these_params, ".csv"))

ggplot(plot_data, aes(x=value, y=variable)) + stat_halfeye() + facet_wrap(~measure, scale="free_x")

pop_dist = read_csv("../BarcelonaTiger/data/proc/biting_networks_veri_census_sections_bcn_pop_props.csv") %>% rename(variable = District, value=vri) %>% mutate(measure="VRI") %>% select(variable, value, measure)

pop_dist_all = pop_dist %>% mutate(variable = "Barcelona (all)")

pop_dist = bind_rows(pop_dist, pop_dist_all)

ggplot(bind_rows(plot_data, pop_dist), aes(x=value, y=variable)) + stat_halfeye() + facet_wrap(~measure, scale="free_x")


hbd = read_csv(paste0("data/sim_summaries/human_bite_distributions", these_params, ".csv")) %>% pivot_longer(cols=everything(), names_to="variable", values_to="value") %>% mutate(measure = "Bites per person") %>% mutate(variable = recode(variable, bcn_all = "Barcelona (all)"))

hbd %>% group_by(variable) %>% summarise(median(value))

H = hbd %>% group_by(variable, k=value) %>% summarise(nk = n()) %>% mutate(N = sum(nk), pk = nk/N, ompksq = (1-pk)^2, dsq = nk*((k-sum(pk))^2)) %>% summarize(N = mean(N), h = sqrt(sum(ompksq)/N), hhet = sqrt(1-3/N), Hm = h/hhet, VAR = mean(sum(dsq))) %>% left_join(plot_data %>% group_by(variable, measure) %>% summarise(mean_value = mean(value)) %>% pivot_wider(id_cols = variable, names_from = measure, values_from = mean_value) )

H %>% ggplot(aes(x = Hm, y = AR)) + geom_point() + geom_smooth(method='lm', formula= y~x)
H %>% ggplot(aes(x = Hm, y = R0)) + geom_point() + geom_smooth(method='lm', formula= y~x)
H %>% ggplot(aes(x = VAR, y = AR)) + geom_point() + geom_smooth(method='lm', formula= y~x)
H %>% ggplot(aes(x = VAR, y = R0)) + geom_point() + geom_smooth(method='lm', formula= y~x)

H %>% dplyr::select(variable, Hm, AR, R0) %>% pivot_longer(cols = -c(variable, Hm), names_to = "outcome") %>% ggplot(aes(x=Hm, y=value)) + geom_point() + geom_smooth(method='lm', formula= y~x) + facet_grid(outcome~., scales = "free_y") + xlab( bquote(Network~heterogeneity~(H[m])))
ggsave(paste0("plots/bcn_districts_hetindex_r0_AR", these_params, ".png"), width=5, height=5)
ggsave(paste0("papers/article_01_biting_networks_scireps/figures/bcn_districts_hetindex_r0_AR", these_params, ".png"), width=5, height=5)


summary(lm(R0~Hm, data = H))

A = bind_rows(lapply(unique(hbd$variable), function(v){
  these_ks = hbd %>% filter(variable == v) %>% pull(value)
  As = 2*sum(abs(dist(these_ks)))
  these_trans_ks = 1/sqrt(these_ks)
  rhos = 2*sum(dist(these_trans_ks))
    return(tibble(variable = v, A = As, rho = rhos))
}))

H %>% left_join(A) %>% ggplot(aes(x = A, y = R0)) + geom_point()  + geom_smooth(method='lm', formula= y~x)
H %>% left_join(A)%>% ggplot(aes(x = A, y = AR)) + geom_point() + geom_smooth(method='lm', formula= y~x)
H %>% left_join(A)%>% ggplot(aes(x = rho, y = R0)) + geom_point() + geom_smooth(method='lm', formula= y~x)
H %>% left_join(A)%>% ggplot(aes(x = rho, y = AR)) + geom_point() + geom_smooth(method='lm', formula= y~x)

H %>% left_join(A) %>% dplyr::select(variable, A, AR, R0) %>% pivot_longer(cols = -c(variable, A), names_to = "outcome") %>% ggplot(aes(x=A, y=value)) + geom_point() + geom_smooth(method='lm', formula= y~x) + facet_grid(outcome~., scales = "free_y") + xlab( bquote(Network~heterogeneity~(A)))
ggsave(paste0("plots/bcn_districts_hetindexA_r0_AR", these_params, ".png"), width=5, height=5)


ggplot(hbd, aes(x=value, y=variable)) + stat_halfeye() #+ facet_wrap(~measure, scale="free_x"))

ggplot(hbd, aes(x=value, color=variable)) + geom_density() #+ facet_wrap(~measure, scale="free_x"))


plot_data = bind_rows(hbd, plot_data, pop_dist) %>% mutate(measure_f = factor(measure, levels=c("VRI", "Bites per person", "R0", "AR")))

these_medians = plot_data %>% group_by(measure_f) %>% summarise(med = median(value))

district_order = plot_data %>% filter(measure == "R0") %>% group_by(variable) %>% summarise(median = median(value)) %>% arrange(median) %>% mutate(order = 1:n())

ggplot(plot_data %>% mutate(variable_f = factor(variable, levels = district_order$variable)), aes(x=value, y=variable_f, fill=measure_f)) + stat_halfeye(normalize="groups", alpha=.7, point_interval = mean_qi,.width = c(0.66, 0.95)) + facet_wrap(~measure_f, scale="free_x", nrow = 1) + theme(legend.position="none") + scale_fill_brewer(palette = "Set1") + xlab("") + ylab("")# + geom_vline(data = these_medians, aes(xintercept=med, color=measure_f)) 
ggsave(paste0("plots/bcn_districts_combo_ordered_vri_hbd_r0_AR", these_params, ".png"), width=8, height=4)


# theoretical distributions ####

combo_data = "data/sim_summaries/combo_plot_AR_R0_eb=800.0_hit=3_mls=20_n_h=100_n_m=400_n_r=2000_n_s=1000_this_set_name=cuel12_tp=0.1.csv"

hd_data = "data/sim_summaries/human_bite_distributions_eb=800.0_hit=3_mls=20_n_h=100_n_m=400_n_r=2000_n_s=1000_this_set_name=cuel12_tp=0.1.csv"

plot_data <- read_csv(combo_data)

ggplot(plot_data, aes(x=value, y=variable)) + stat_halfeye() + facet_wrap(~measure, scale="free_x")


hbd = read_csv(hd_data) %>% pivot_longer(cols=everything(), names_to="variable", values_to="value") %>% mutate(measure = "Bites per person")

# THOUGHT: ADD 1 to k to avoid zeros?

H = hbd %>% group_by(variable, k=value) %>% summarise(nk = n()) %>% mutate(N = sum(nk), pk = nk/N, ompksq = (1-pk)^2, dsq = nk*((k-sum(pk))^2)) %>% summarize(N = mean(N), h = sqrt(sum(ompksq)/N), hhet = sqrt(1-3/N), Hm = h/hhet, VAR = mean(sum(dsq))) %>% left_join(plot_data %>% group_by(variable, measure) %>% summarise(mean_value = mean(value)) %>% pivot_wider(id_cols = variable, names_from = measure, values_from = mean_value) )

H %>% ggplot(aes(x = Hm, y = AR)) + geom_point()
H %>% ggplot(aes(x = Hm, y = R0)) + geom_point()
H %>% ggplot(aes(x = VAR, y = AR)) + geom_point()
H %>% ggplot(aes(x = VAR, y = R0)) + geom_point()

summary(lm(AR~Hm, data = H))

summary(lm(R0~Hm, data = H))


ggplot(hbd, aes(x=value, y=variable)) + stat_halfeye() #+ facet_wrap(~measure, scale="free_x"))

ggplot(hbd, aes(x=value, color=variable)) + geom_density() #+ facet_wrap(~measure, scale="free_x"))


unique(plot_data$variable)

plot_data = bind_rows(hbd, plot_data) %>% mutate(measure_f = factor(measure, levels=c("Bites per person", "R0", "AR")), variable_f = factor(variable, levels=c("exp", "uniform", "levy1p1", "levy2p1", "constant"), labels = c("exponenaial", "uniform", "levy 1.1", "levy 2.1", "constant")))

these_medians = plot_data %>% group_by(measure_f) %>% summarise(med = median(value))

this_pal = brewer.pal(4, "Set1")[2:4]

ggplot(plot_data, aes(x=value, y=variable_f, fill=measure_f)) + stat_halfeye(normalize="groups", alpha=.7, point_interval = mean_qi,.width = c(0.66, 0.95)) + facet_wrap(~measure_f, scale="free_x", nrow = 1) + theme(legend.position="none") + scale_fill_manual(values=this_pal) + xlab("") + ylab("")# + geom_vline(data = these_medians, aes(xintercept=med, color=measure_f)) 
# ggsave("plots/theoretical_dist_combo_vri_hbd_r0_AR.png", width=8, height=4)


A = bind_rows(lapply(unique(hbd$variable), function(v){
  these_ks = hbd %>% filter(variable == v) %>% pull(value)
  As = 2*sum(abs(dist(these_ks)))
  these_trans_ks = 1/sqrt(these_ks+1) # NOTE I AM ADDING 1 HERE TO AVOID ZEROS IN DENOMINATOR. THINK IF THIS IS OK...
  rhos = 2*sum(dist(these_trans_ks))
  return(tibble(variable = v, A = As, rho = rhos))
}))

H %>% left_join(A) %>% ggplot(aes(x = A, y = R0)) + geom_point()  + geom_smooth(method='lm', formula= y~x)
H %>% left_join(A)%>% ggplot(aes(x = A, y = AR)) + geom_point() + geom_smooth(method='lm', formula= y~x)
H %>% left_join(A)%>% ggplot(aes(x = rho, y = R0)) + geom_point() + geom_smooth(method='lm', formula= y~x)
H %>% left_join(A)%>% ggplot(aes(x = rho, y = AR)) + geom_point() + geom_smooth(method='lm', formula= y~x)

H %>% left_join(A) %>% dplyr::select(variable, A, AR, R0) %>% pivot_longer(cols = -c(variable, A), names_to = "outcome") %>% ggplot(aes(x=A, y=value)) + geom_point() + geom_smooth(method='lm', formula= y~x) + facet_grid(outcome~., scales = "free_y") + xlab( bquote(Network~heterogeneity~(A)))
#ggsave(paste0("plots/bcn_districts_hetindexA_r0_AR", these_params, ".png"), width=5, height=5)


