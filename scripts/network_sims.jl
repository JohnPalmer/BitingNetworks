using Base: Multimedia
using Random, Distributions

Random.seed!(123)

latent_period = 3
recovery_period = 3

transmission_prob = 1

n_humans = 50
mean_bites_per_person = 10

n_mosquitoes = 100
mean_bites_per_mosquito = 2

human_distribution = Beta(.5, .5)
mosquito_distribution = Beta(.5, .5)

human_probs = rand(human_distribution, n_humans)
mosquito_probs = rand(mosquito_distribution, n_mosquitoes)

# creating first infection
status_humans = zeros(Int8, n_humans)
status_humans[rand(1:n_humans, 1)[1]] = 1

status_mosquitoes = zeros(Int8, n_mosquitoes)

# human-to-mosquito infections
for h in findall(status_humans .== 1), m in findall(status_mosquitoes .==0)
  if rand(Bernoulli(human_probs[h]*mosquito_probs[m]*transmission_prob), 1)[1]
    status_mosquitoes[m] = 1
  end
end

# mosquito-human infections
for h in findall(status_humans .== 0), m in findall(status_mosquitoes .==1)
  if rand(Bernoulli(human_probs[h]*mosquito_probs[m]*transmission_prob), 1)[1]
    status_humans[h] = 1
  end
end


