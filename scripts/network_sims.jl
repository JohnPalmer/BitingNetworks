using Random, Distributions
using ProgressMeter, Plots
gr()

Threads.nthreads()

Random.seed!(123)

n_steps = 200
n_humans = 500
mean_bites_per_person = 10
n_mosquitoes = 1000
mean_bites_per_mosquito = 2

latent_period = 3
recovery_period = 3

transmission_prob = .0005

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

  i_hs = findall(status_humans .== 1)
  s_ms = findall(status_mosquitoes .==0)
  these_iterations = 1:length(i_hs)*length(s_ms)

# human-to-mosquito infections
  Threads.@threads for i in these_iterations
    h = i_hs[cld(i, length(s_ms))]
    m = s_ms[mod1(i, length(s_ms))]
    if rand(Bernoulli(human_probs[h]*mosquito_probs[m]*transmission_prob), 1)[1]
      status_mosquitoes[m] = 1
    end
  end

n_mosquito_infections[s] = sum(status_mosquitoes)

s_hs = findall(status_humans .== 0)
i_ms = findall(status_mosquitoes .==1)
these_iterations = 1:length(s_hs)*length(i_ms)

  # mosquito-human infections
  Threads.@threads for i in these_iterations
    h = s_hs[cld(i, length(i_ms))]
    m = i_ms[mod1(i, length(i_ms))]
    if rand(Bernoulli(human_probs[h]*mosquito_probs[m]*transmission_prob), 1)[1]
      status_humans[h] = 1
    end
  end

  n_human_infections[s] = sum(status_humans)

end

using Plots
gr()


infection_data = [n_human_infections, n_mosquito_infections]

infection_data_props = [n_human_infections/n_humans, n_mosquito_infections/n_mosquitoes]

plot(1:n_steps, infection_data_props, lab=["humans" "mosquitoes"], w=3, legend=:right)


