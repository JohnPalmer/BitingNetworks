using Random, Distributions
using ProgressMeter, Plots, StatsPlots, KernelDensity, Bootstrap
gr()

using Bites

Threads.nthreads()

Random.seed!(123)

n_steps = 50
n_reps = 300
n_humans = 1000
n_mosquitoes = 4000

transmission_prob = .4

expected_bites = n_mosquitoes*2.0

# how long before infected human recovers
const human_infection_time = 3

# how long do mosquitoes live. Since mosquitoes do not usually get accute arbovirus infections, the life span is the key: They will go from susceptible->infected->dead/reborn. (With assumption of population stability, deaths simply put the agent back to susceptible status.) "Unlike arboviral infections in humans, which are usually acute, arboviral infections in mosquitoes are persistent. Once the infection is established, the mosquito remains infected for the rest of its life." https://www.sciencedirect.com/science/article/pii/S1931312819303701
const mosquito_life_span = 20

HD_const, MD_const = distribute_bite_probabilities( (1/n_mosquitoes):(1/n_mosquitoes), (1/n_humans):(1/n_humans), n_humans, n_mosquitoes, expected_bites)


HD_uniform_0_1, MD_const = distribute_bite_probabilities(Distributions.Uniform(0,1), (1/n_humans):(1/n_humans), n_humans, n_mosquitoes, expected_bites)

HD_levy_1_0001, MD_const = distribute_bite_probabilities(Distributions.Levy(1, .0001), (1/n_humans):(1/n_humans), n_humans, n_mosquitoes, expected_bites)

HD_exp05, MD_const = distribute_bite_probabilities(Distributions.Exponential(1/.5), (1/n_humans):(1/n_humans), n_humans, n_mosquitoes, expected_bites)


HD_truncnormal_1_3 = distribute_bite_probabilities(Truncated(Normal(1, 3), 0, Inf), (1/n_humans):(1/n_humans), n_humans, n_mosquitoes, expected_bites)[1]


HD_trunclevy_1_3 = distribute_bite_probabilities(Truncated(Levy(1.1, 3), 1, 1000), (1/n_humans):(1/n_humans), n_humans, n_mosquitoes, expected_bites)[1]


scenario_results = []

scenarios = ( 
  constant= [MD_const, HD_const],
  uniform = [MD_const, HD_uniform_0_1],
  trunclevy = [MD_const, HD_trunclevy_1_3],
  exp = [MD_const, HD_exp05],
  levy = [MD_const, HD_levy_1_0001])


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

human_R0_converg_check = Vector(undef, n_reps)
for i in 1:n_reps
  human_R0_converg_check[i] = mean(R0s[1:i])
end

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

mosquito_R0_converg_check = Vector(undef, n_reps)
for i in 1:n_reps
  mosquito_R0_converg_check[i] = mean(R0s[1:i])
end


scenario_result = (n_human_infections_reps=n_human_infections_reps, human_R0=human_R0, mosquito_R0=mosquito_R0, human_bite_distribution=human_bite_distribution, mosquito_bite_distribution=mosquito_bite_distribution, n_human_recovered_reps=n_human_recovered_reps, human_R0_converg_check=human_R0_converg_check, mosquito_R0_converg_check=mosquito_R0_converg_check)

push!(scenario_results, scenario_result)

end

human_R0_convergence_checks = [mean(x[:human_R0_converg_check][1:i]) for x in scenario_results, i in 1:n_reps]

R0 = [x[2] for x in scenario_results]

n_boot = 1000
human_R0_boot = [stderror(bootstrap(mean, x[:human_R0_converg_check], BasicSampling(n_boot))) for x in scenario_results]

these_labs = [string(x) for x in keys(scenarios)]

plot(transpose(human_R0_convergence_checks), label=reshape(these_labs, 1, length(scenario_results)), xlabel="Repititions", ylabel="Mean R0", legend=:bottomright, linewidth=2)
plot!(size=(800,600))
png("plots/R0_convergence_check_plot.png")


AR = [ [mean(x[6][:,n_steps] .+ x[1][:,n_steps]), quantile(x[6][:,n_steps] .+ x[1][:,n_steps], [.025, .25, .5, .75, .975])] for x in scenario_results]

AR_convergence_check = [mean(x[6][1:i,n_steps] .+ x[1][1:i,n_steps]) for x in scenario_results, i in 1:n_reps]

plot(transpose(AR_convergence_check), label=reshape([string("S", x) for x in 1:length(scenario_results)], 1, length(scenario_results)), xlabel="Repititions", ylabel="Mean Attack Rate")
plot!(size=(800,600))
png("plots/AR_convergence_check_plot.png")



mean_bites_per_person = [mean(x[4]) for x in scenario_results]

max_bites_per_person = [maximum(x[4]) for x in scenario_results]


maximum(scenario_results[1][1])

these_plots = []

this_max = maximum([maximum(x[1]) for x in scenario_results])

steps_to_plot = 15

for i in 1:length(scenarios)
  n_human_infections_reps = scenario_results[i][1]./n_humans
  n_human_infections_reps_median = transpose(median(n_human_infections_reps, dims=1))

  this_p = plot(1:steps_to_plot, transpose(n_human_infections_reps[ :, 1:steps_to_plot]), w=.2, color=:grey, legend=:none, yaxis = ("Infected", (0, 1), 0:.5:1), xaxis = ("Time"))

  plot!(this_p, 1:steps_to_plot,n_human_infections_reps_median[1:steps_to_plot], w=2, color=:black, yaxis = ("Infected", (0, 1), 0:.2:1))

  push!(these_plots, this_p)
end

max_bites = maximum([maximum(x[4]) for x in scenario_results])

kd_max = maximum(kde(vec(scenario_results[1][4])).density)

these_bite_plots = [density(scenario_results[i][:human_bite_distribution], xaxis = ("N Bites", (0, max_bites), 0:50:max_bites), yaxis = ("Density", (-0.01, kd_max)), legend=:none, normed=true, linewidth=2, title=these_labs[i]) for i in 1:length(scenario_results)]

l = @layout [a f; b g; c h; d i; e j]

plot(these_plots[1], these_bite_plots[1], these_plots[2], these_bite_plots[2], these_plots[3], these_bite_plots[3],these_plots[4], these_bite_plots[4], these_plots[5], these_bite_plots[5], layout = l)
plot!(size=(800,800))
png("plots/epicurve_comparison.png")