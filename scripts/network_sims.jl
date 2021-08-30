using Random, Distributions
using ProgressMeter, Plots, StatsPlots, KernelDensity, Bootstrap
plotlyjs()

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

fixed_mosquito_probs = (1/n_humans):(1/n_humans)

subpop_a_size = Int(n_humans/2)
subpop_a = rand(Truncated(Normal(3, 3), 0, Inf), subpop_a_size)

subpop_b_size = n_humans - subpop_a_size
subpop_b = rand(Truncated(Normal(30, 3), 0, Inf), subpop_b_size)

rand(subpop_b, 10)

scenario_results = []

human_distributions = ( 
  constant= (1/n_mosquitoes):(1/n_mosquitoes),
  uniform = Uniform(0,1),
  exp = Exponential(1/.5),
  tlevy1p1 = Truncated(Levy(1.1, .0001), 1, 1000),
  levy1p1 = Levy(1.1, .0001),
  levy2p1 = Levy(2.1, .0001),
  levy3p1 = Levy(3.1, .0001),
  mixed_norms = vcat(subpop_a, subpop_b)
)


  human_distribution = human_distributions[8]

for human_distribution in human_distributions

  human_probs, mosquito_probs = distribute_bite_probabilities(human_distribution, fixed_mosquito_probs, n_humans, n_mosquitoes, expected_bites)

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

# density(human_bite_distribution)

  n_human_infections_reps = Array{Int64}(undef, n_reps, n_steps)

  n_mosquito_infections_reps = Array{Int64}(undef, n_reps, n_steps)

  n_human_recovered_reps = Array{Int64}(undef, n_reps, n_steps)

  r = 1

  p = Progress(n_reps)

  Threads.@threads for r = 1:n_reps

    human_probs, mosquito_probs = distribute_bite_probabilities(human_distribution, fixed_mosquito_probs, n_humans, n_mosquitoes, expected_bites)
  
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
    if sum(x.>this_min)>0
      this_next_min = minimum(x[x.>this_min])
      R0s[i] = (this_next_min - this_min) / this_min
    else
      R0s[i] = 0
    end
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

these_labs = [string(x) for x in keys(human_distributions)]

plot(transpose(human_R0_convergence_checks), label=reshape(these_labs, 1, length(scenario_results)), xlabel="Repetitions", ylabel="Mean R0", legend=:bottomright, linewidth=2, palette = :Dark2_8)
plot!(size=(800,600))
png("plots/R0_convergence_check_plot.png")


AR = [ [mean(x[6][:,n_steps] .+ x[1][:,n_steps]), quantile(x[6][:,n_steps] .+ x[1][:,n_steps], [.025, .25, .5, .75, .975])] for x in scenario_results]

AR_convergence_check = [mean(x[6][1:i,n_steps] .+ x[1][1:i,n_steps]) for x in scenario_results, i in 1:n_reps]

plot(transpose(AR_convergence_check), label=label=reshape(these_labs, 1, length(scenario_results)), xlabel="Repititions", ylabel="Mean Attack Rate", palette = :Dark2_8, linewidth=2)
plot!(size=(800,600))
png("plots/AR_convergence_check_plot.png")



mean_bites_per_person = [mean(x[4]) for x in scenario_results]

max_bites_per_person = [maximum(x[4]) for x in scenario_results]


maximum(scenario_results[1][1])

using Plots.PlotMeasures
these_plots = []

this_max = maximum([maximum(x[1]) for x in scenario_results])

steps_to_plot = 15

for i in 1:length(human_distributions)
  n_human_infections_reps = scenario_results[i][1]./n_humans
  n_human_infections_reps_median = transpose(median(n_human_infections_reps, dims=1))

  this_p = plot(1:steps_to_plot, transpose(n_human_infections_reps[ :, 1:steps_to_plot]), w=.2, color=:grey, legend=:none, yaxis = ("Infected", (0, 1), 0:.5:1), xaxis = ("Time"), margin = 1cm)


  plot!(this_p, 1:steps_to_plot,n_human_infections_reps_median[1:steps_to_plot], w=2, color=:black, yaxis = ("Infected", (0, 1), 0:.2:1))

  push!(these_plots, this_p)
end

max_bites = maximum([maximum(x[4]) for x in scenario_results])

kd_max = maximum(kde(vec(scenario_results[1][4])).density)

these_bite_plots = [density(scenario_results[i][:human_bite_distribution], xaxis = ("N Bites", (0, max_bites), 0:ceil(max_bites/10):max_bites), yaxis = ("Density", (-0.01, kd_max)), legend=:none, normed=true, linewidth=2, title=these_labs[i], margin = 1cm) for i in 1:length(scenario_results)]


this_p = plot( 
  these_bite_plots[1], 
  these_plots[1], 
  these_bite_plots[2], 
  these_plots[2], 
  these_bite_plots[3], 
  these_plots[3], 
  these_bite_plots[4], 
  these_plots[4], 
  these_bite_plots[5], 
  these_plots[5], 
  these_bite_plots[6], 
  these_plots[6], 
  these_bite_plots[7], 
  these_plots[7], 
  these_bite_plots[8], 
  these_plots[8], 
  layout=(8,2)
)
plot!(this_p, size=(800,2000))
savefig(this_p, "plots/epicurve_comparison.html")