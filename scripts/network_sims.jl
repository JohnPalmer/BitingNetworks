using Random, Distributions
using ProgressMeter, Plots
gr()

Random.seed!(123)

n_steps = 100
n_humans = 5000
mean_bites_per_person = 10
n_mosquitoes = 10000
mean_bites_per_mosquito = 2

latent_period = 3
recovery_period = 3

transmission_prob = .001

human_distribution = Beta(.5, .5)
mosquito_distribution = Beta(.5, .5)

human_probs = rand(human_distribution, n_humans)
mosquito_probs = rand(mosquito_distribution, n_mosquitoes)

# creating first infection
status_humans = zeros(Int8, n_humans)
status_humans[rand(1:n_humans, 1)[1]] = 1

status_mosquitoes = zeros(Int8, n_mosquitoes)

n_human_infections = Vector{Int64}(undef, n_steps)
n_mosquito_infections = Vector{Int64}(undef, n_steps)

n_human_infections[1] = 1
n_mosquito_infections[1] = 0

@showprogress for s = 2:n_steps

# human-to-mosquito infections
  for h in findall(status_humans .== 1), m in findall(status_mosquitoes .==0)
    if rand(Bernoulli(human_probs[h]*mosquito_probs[m]*transmission_prob), 1)[1]
      status_mosquitoes[m] = 1
    end
  end

n_mosquito_infections[s] = sum(status_mosquitoes)

  # mosquito-human infections
  for h in findall(status_humans .== 0), m in findall(status_mosquitoes .==1)
    if rand(Bernoulli(human_probs[h]*mosquito_probs[m]*transmission_prob), 1)[1]
      status_humans[h] = 1
    end
  end

  n_human_infections[s] = sum(status_humans)

end

using Plots
gr()

plot(1:n_steps, n_human_infections)

