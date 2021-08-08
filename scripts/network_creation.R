library(tidyverse)
library(network)
library(bipartite) 

n_humans = 10
mean_bites_per_person = 10

n_mosquitoes = 100
mean_bites_per_mosquito = 2

human_probs = rpois(n_humans, mean_bites_per_person)/n_mosquitoes

mosquito_probs =rpois(n_mosquitoes, mean_bites_per_mosquito)/n_humans

bn_tab = human_probs %*% t(mosquito_probs) %>% as_tibble(rownames = "human_id") %>% pivot_longer(cols = -human_id, names_to = "mosquito_id", values_to = "bite_prob") 

bn_tab$bite = rbinom(nrow(bn_tab), size = 1, bn_tab$bite_prob)
  
bn_tab %>% group_by(human_id) %>% summarise(sum(bite))

bn_tab %>% group_by(mosquito_id) %>% summarise(sum(bite))
