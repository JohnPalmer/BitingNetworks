library(tidyverse)
library(network)
library(bipartite) 
library(data.table)

latent_period = 3
recovery_period = 3

n_humans = 5
mean_bites_per_person = 10
human_ids = paste("H", 1:n_humans)

n_mosquitoes = 10
mean_bites_per_mosquito = 2
mosquito_ids = paste("M", 1:n_mosquitoes)

human_probs = rpois(n_humans, mean_bites_per_person)/n_mosquitoes

mosquito_probs =rpois(n_mosquitoes, mean_bites_per_mosquito)/n_humans

# human_probs = runif(n_humans, min=0, max = mean_bites_per_person*2)/n_mosquitoes

# mosquito_probs = runif(n_mosquitoes, min=0, max = mean_bites_per_mosquito*2)/n_humans

# nbinom_dispersion = 10

# human_probs = rnbinom(n=n_humans, size = nbinom_dispersion, mu= mean_bites_per_person)/n_mosquitoes

# mosquito_probs =rnbinom(n = n_mosquitoes,  size = nbinom_dispersion, mu = mean_bites_per_mosquito)/n_humans

# human_probs = rep(mean_bites_per_person, n_humans)/n_mosquitoes

# mosquito_probs = rep(mean_bites_per_mosquito, n_mosquitoes)/n_humans


human_infections = tibble(human_id = human_ids, infected_human = c(1, rep(0, n_humans-1)))

mosquito_infections = tibble(mosquito_id = mosquito_ids, infected_mosquito = 0)

# matrix of probabilities of each mosquito biting each person. People are the rows, mosquitoes are the columns
prob_mat = human_probs %*% t(mosquito_probs) 

# For each human-mosquito pair I then take a draw from the binomial distribution to see if there is a bite
bite_mat = matrix(rbinom(n = length(prob_mat), size = 1, prob = as.vector(prob_mat)), nrow = dim(prob_mat)[1], ncol=dim(prob_mat)[2])

# creating infections
infected_humans = rep(0, n_humans)
infected_humans[sample(1:n_humans, 1)] = 1
infected_mosquitoes = rep(0, n_mosquitoes)

new_infections = infected_humans*bite_mat + infected_mosquitoes*bite_mat

infected_humans[which(rowSums(new_infections) > 0)] = 1

infected_mosquitoes[which(colSums(new_infections) > 0)] = 1

infection_mat * bite_mat


#%>% as_tibble("human_id" = 1:n_humans) %>% pivot_longer(cols = -human_id, names_to = "mosquito_id", values_to = "bite_prob") %>% left_join(human_infections) %>% left_join(mosquito_infections) %>% as.data.table()


bn_tab$bite = rbinom(nrow(bn_tab), size = 1, bn_tab$bite_prob)
  
bn_tab %>% group_by(human_id) %>% summarise(sum(bite))

bn_tab %>% group_by(mosquito_id) %>% summarise(sum(bite))


# dynamics ####

n_steps = 100

infection_summary = tibble(step = 1:n_steps, humans = rep(0, n_steps), mosquitoes = rep(0, n_steps))

i = 1

for(i in 1:n_steps){
  
  bn_tab[bite==1, infected_human:=as.integer((infected_human+infected_mosquito)>0)]
  bn_tab[bite==1, infected_mosquito:=as.integer((infected_human+infected_mosquito)>0)]
  
  human_infections = bn_tab %>% group_by(human_id) %>% summarize(infected_human = max(infected_human))
  
  mosquito_infections = bn_tab %>% group_by(mosquito_id) %>% summarize(infected_mosquito = max(infected_mosquito))
  
  infection_summary$humans[i] = sum(human_infections$infected_human)
  
  infection_summary$mosquitoes[i] = sum(mosquito_infections$infected_mosquito)

  # update
  
  bn_tab = human_probs %*% t(mosquito_probs) %>% as_tibble(rownames = "human_id") %>% pivot_longer(cols = -human_id, names_to = "mosquito_id", values_to = "bite_prob") %>% left_join(human_infections) %>% left_join(mosquito_infections) %>% as.data.table()

  
  bn_tab$bite = rbinom(nrow(bn_tab), size = 1, bn_tab$bite_prob)
  

}

infection_summary_long = infection_summary %>% pivot_longer(cols = -step, names_to = "species", values_to = "value") %>% as.data.table()

infection_summary_long$percent = 0
infection_summary_long[species=="mosquitoes", percent:=value/n_mosquitoes]
infection_summary_long[species=="humans", percent:=value/n_humans]

ggplot(infection_summary_long, aes(x=step, y=percent, color=species)) + geom_line()


daily_summary = tibble(step = 1:(n_steps-1), humans = diff(infection_summary$humans), mosquitoes = diff(infection_summary$mosquitoes))

daily_summary_long = daily_summary %>% pivot_longer(cols = -step, names_to = "species", values_to = "value") %>% as.data.table()


ggplot(daily_summary_long, aes(x=step, y=value, color=species)) + geom_line()
