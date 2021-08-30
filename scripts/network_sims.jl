using Random, Distributions
using ProgressMeter, Plots, StatsPlots
gr()

using Bites

Threads.nthreads()

Random.seed!(123)

n_steps = 50
n_reps = 300
n_humans = 1000
n_mosquitoes = 4000

transmission_prob = .06

# how long before infected human recovers
const human_infection_time = 3

# how long do mosquitoes live. Since mosquitoes do not usually get accute arbovirus infections, the life span is the key: They will go from susceptible->infected->dead/reborn. (With assumption of population stability, deaths simply put the agent back to susceptible status.) "Unlike arboviral infections in humans, which are usually acute, arboviral infections in mosquitoes are persistent. Once the infection is established, the mosquito remains infected for the rest of its life." https://www.sciencedirect.com/science/article/pii/S1931312819303701
const mosquito_life_span = 20

const neg_bin_r = 2
const neg_bin_p = .1

const neg_bin_r_a = 5
const neg_bin_p_a = .21

const neg_bin_r_b = 10
const neg_bin_p_b = .35

human_distribution = NegativeBinomial(neg_bin_r, neg_bin_p)

human_distribution_a = NegativeBinomial(neg_bin_r_a, neg_bin_p_a)

human_distribution_b = NegativeBinomial(neg_bin_r_b, neg_bin_p_b)

# human_distribution = Beta(2, 5)
# mosquito_distribution = Beta(.2, 2)
mosquito_probs1 = fill(1/n_humans, n_mosquitoes)

human_probs1 = rand(human_distribution, n_humans)/(sum(mosquito_probs1))

human_probs1a = rand(human_distribution_a, n_humans)/(sum(mosquito_probs1))

human_probs1b = rand(human_distribution_b, n_humans)/(sum(mosquito_probs1))

human_probs2 = fill(mean(human_probs1), n_humans)
mosquito_probs2 = fill(mean(mosquito_probs1), n_mosquitoes)


human_distribution3 = Uniform(0, mean(human_probs1)*2)
mosquito_distribution3 = Uniform(0, mean(mosquito_probs1)*2)

human_probs3 = rand(human_distribution3, n_humans)
mosquito_probs3 = rand(mosquito_distribution3, n_mosquitoes)


scenario_results = []

scenarios = [[mosquito_probs1, human_probs1],[mosquito_probs1, human_probs1a], [mosquito_probs1, human_probs1b],[mosquito_probs2, human_probs2], [mosquito_probs3, human_probs3] ]

scenario = scenarios[4]

for scenario in scenarios

  mosquito_probs = scenario[1]

  human_probs = scenario[2]

  these_iterations = 1:length(human_probs)*length(mosquito_probs)

  bn = BitArray(undef, length(human_probs), length(mosquito_probs))
  # biting network snapshot
  Threads.@threads for i in these_iterations
    this_row = cld(i, length(mosquito_probs))
    this_col = mod1(i, length(mosquito_probs))
    h = human_probs[this_row]
    m = mosquito_probs[this_col]
    joint_prob = h*m
    bn[this_row, this_col] = rand(Bernoulli( ifelse(joint_prob<1, joint_prob, 1)), 1)[1]
  end


  mosquito_bite_distribution = transpose(sum(bn, dims=1))

  human_bite_distribution = sum(bn, dims=2)

#  density(human_bite_distribution)

#  mean(human_bite_distribution)

#  density(mosquito_bite_distribution)

  n_human_infections_reps = Array{Int64}(undef, n_reps, n_steps)

  n_mosquito_infections_reps = Array{Int64}(undef, n_reps, n_steps)

  n_human_recovered_reps = Array{Int64}(undef, n_reps, n_steps)

  r = 1

  p = Progress(n_reps)

  Threads.@threads for r = 1:n_reps


    n_mosquito_infections_reps[r, :], n_human_infections_reps[r, :], n_human_recovered_reps[r, :] = bite_steps(n_steps, n_humans, n_mosquitoes, human_infection_time, mosquito_life_span, human_probs, mosquito_probs, transmission_prob)


    next!(p)

  end

  R0s = Vector(undef, n_reps)

  for i in 1:n_reps
    x = n_human_infections_reps[i,:]
    if maximum(x) > 1
      R0s[i] = minimum(x[x.>1])
    else
      R0s[i] = 0
    end
  end

human_R0 = mean(R0s)


R0s = Vector(undef, n_reps)

for i in 1:n_reps
  x = n_mosquito_infections_reps[i,:]
  if maximum(x) > 1
    this_min = minimum(x[x.>0])
    this_next_min = minimum(x[x.>this_min])
    R0s[i] = (this_next_min - this_min) / this_min
  else
    R0s[i] = 0
  end
end

mosquito_R0 = mean(R0s)

scenario_result = (n_human_infections_reps=n_human_infections_reps, human_R0=human_R0, mosquito_R0=mosquito_R0, human_bite_distribution=human_bite_distribution, mosquito_bite_distribution=mosquito_bite_distribution, n_human_recovered_reps=n_human_recovered_reps)

push!(scenario_results, scenario_result)

end

R0 = [x[2] for x in scenario_results]

AR = [ [mean(x[6][:,n_steps] .+ x[1][:,n_steps]), quantile(x[6][:,n_steps] .+ x[1][:,n_steps], [.025, .25, .5, .75, .975])]  for x in scenario_results]

mean_bites_per_person = [mean(x[4]) for x in scenario_results]


these_plots = []

this_max = maximum([maximum(x[1]) for x in scenario_results])

this_max = 100*cld(this_max, 100)

for i in 1:length(scenarios)
  n_human_infections_reps = scenario_results[i][1]
  n_human_infections_reps_median = transpose(median(n_human_infections_reps, dims=1))

  this_p = plot(1:n_steps, transpose(n_human_infections_reps), w=.2, color=:grey, legend=:none, yaxis = ("Infected", (0, this_max), 0:100:this_max), xaxis = ("Time"))

  plot!(this_p, 1:n_steps,n_human_infections_reps_median, w=2, color=:black, yaxis = ("Infected", (0, this_max), 0:100:this_max))

  push!(these_plots, this_p)
end

max_bites = maximum([maximum(x[4]) for x in scenario_results])

these_bite_plots = [density(x[4], xaxis = ("N Bites", (0, max_bites), 0:40:max_bites), yaxis = ("Density", (0, .09)), legend=:none) for x in scenario_results]

l = @layout [a f; b g; c h; d i; e j]

plot(these_plots[1], these_bite_plots[1], these_plots[2], these_bite_plots[2], these_plots[3], these_bite_plots[3],these_plots[4], these_bite_plots[4], these_plots[5], these_bite_plots[5], layout = l)
