library(tidyverse)
library(jsonlite)
library(lubridate)

reps = read_csv("~/research/data/mosquito_alert_data/tigaserver_app_report.csv")
reps %>% filter(type == "bite") %>% pull(creation_time) %>% range()

D = read_csv("~/research/data/mosquito_alert_data/tigaserver_app_reportresponse.csv")

n_bites = D %>% filter(!is.na(answer_value) & question_id==1)

mean(n_bites$answer_value)
sd(n_bites$answer_value)
range(n_bites$answer_value)
mean(n_bites$answer_value>3)
n_bites %>% ggplot(aes(x=answer_value)) + geom_density()

ggplot(n_bites, aes(x=answer_value)) +
  stat_density(aes(y=..count..), color="black", fill="blue", alpha=0.3) +
  scale_x_continuous(trans="log1p", expand=c(0,0)) + theme_bw() + xlab("Number of bites") + ylab("Number of people")

