using DrWatson
quickactivate(@__DIR__, "BitingNetworks")
using Bites, Random, Distributions
using ProgressMeter, Plots, StatsPlots, KernelDensity, Bootstrap
plotlyjs()

Threads.nthreads()

Random.seed!(123)

n_steps = 50
n_reps = 600
n_humans = 10000
n_mosquitoes = 40000

transmission_prob = .2

expected_bites = float(n_mosquitoes)

# how long before infected human recovers
const human_infection_time = 3

# how long do mosquitoes live. Since mosquitoes do not usually get accute arbovirus infections, the life span is the key: They will go from susceptible->infected->dead/reborn. (With assumption of population stability, deaths simply put the agent back to susceptible status.) "Unlike arboviral infections in humans, which are usually acute, arboviral infections in mosquitoes are persistent. Once the infection is established, the mosquito remains infected for the rest of its life." https://www.sciencedirect.com/science/article/pii/S1931312819303701
const mosquito_life_span = 20

fixed_mosquito_prob = (1/n_humans)
mosquito_distribution = fixed_mosquito_prob:fixed_mosquito_prob
mosquito_distribution_name = string("fixed_", fixed_mosquito_prob)

subpop_a_size = Int(n_humans/2)
subpop_a_distribution = Truncated(Normal(3, 3), 0, Inf)
subpop_a = rand(subpop_a_distribution, subpop_a_size)

subpop_b_size = n_humans - subpop_a_size
subpop_b_distribution = Truncated(Normal(30, 3), 0, Inf)
subpop_b = rand(subpop_b_distribution, subpop_b_size)

scenario_results = []

human_distributions = ( 
  constant = (1/n_mosquitoes):(1/n_mosquitoes),
  uniform = Uniform(0,1),
  exp = Exponential(1/.5),
  tlevy1p1 = Truncated(Levy(1.1, .0001), 1, 1000),
  tlevy2p1 = Truncated(Levy(2.1, .0001), 1, 1000),
  tlevy3p1 = Truncated(Levy(3.1, .0001), 1, 1000),
  mixed_norms = vcat(subpop_a, subpop_b)
)

human_distribution_names = join([string(x) for x in keys(human_distributions)], ".")

this_sim_dict_a = @strdict n_steps n_reps n_humans n_mosquitoes transmission_prob expected_bites human_infection_time mosquito_life_span

this_sim_dict = @strdict n_steps n_reps n_humans n_mosquitoes transmission_prob expected_bites human_infection_time mosquito_life_span human_distribution_names

hdi = 1

for hdi in 1:length(human_distributions)

  human_distribution = human_distributions[hdi]
  human_distribution_name = string(keys(human_distributions)[hdi])

  human_probs, mosquito_probs = distribute_bite_probabilities(human_distribution, mosquito_distribution, n_humans, n_mosquitoes, expected_bites)

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

    human_probs, mosquito_probs = distribute_bite_probabilities(human_distribution, mosquito_distribution, n_humans, n_mosquitoes, expected_bites)
  
    n_mosquito_infections_reps[r, :], n_human_infections_reps[r, :], n_human_recovered_reps[r, :] = bite_steps(n_steps, n_humans, n_mosquitoes, human_infection_time, mosquito_life_span, human_probs, mosquito_probs, transmission_prob)

    next!(p)

  end

  human_R0_results = calculate_r0(n_reps, n_human_infections_reps, 1)
  mosquito_R0_results = calculate_r0(n_reps, n_mosquito_infections_reps, 0)

scenario_result = (
  n_human_infections_reps = n_human_infections_reps, 
  human_R0=human_R0_results.R0, 
  mosquito_R0 = mosquito_R0_results.R0, 
  human_bite_distribution = human_bite_distribution, 
  mosquito_bite_distribution = mosquito_bite_distribution,
  n_human_recovered_reps = n_human_recovered_reps,
  human_R0_converge_check = human_R0_results.converge_check,
  mosquito_R0_converge_check = mosquito_R0_results.converge_check,
  human_R0_reps = human_R0_results.R0_reps,
  mosquito_R0_reps= mosquito_R0_results.R0_reps, name = human_distribution_name,
  n_steps = 50,
  n_reps = 300,
  n_humans = 1000,
  n_mosquitoes = 4000,
  trans_prob = .4, 
  exp_bites = n_mosquitoes*2.0,
  human_infect_time = 3,
  mosq_life_span = 20,
  human_dist = human_distribution,
  human_dist_name = human_distribution_name,
  mosquito_dist = mosquito_distribution,
  mosquito_dist_name = mosquito_distribution_name
  )

@tagsave(datadir("sims", savename(scenario_result, "jld2")), tostringdict(ntuple2dict(scenario_result)))

push!(scenario_results, scenario_result)

end


####
human_R0_convergence_checks = [x[:human_R0_converge_check] for x in scenario_results]


R0 = [x.human_R0 for x in scenario_results]

these_labs = [string(x) for x in keys(human_distributions)]

plot(human_R0_convergence_checks, label=reshape(these_labs, 1, length(scenario_results)), xlabel="Repetitions", ylabel="Mean R0", legend=:bottomright, linewidth=2, palette = :Dark2_8)
plot!(size=(800,600))
png(string("plots/R0_conv", savename(this_sim_dict), ".png"))


AR = [ [mean(x[:n_human_recovered_reps][:,n_steps] .+ x[:n_human_infections_reps][:,n_steps]), quantile(x[6][:,n_steps] .+ x[1][:,n_steps], [.025, .25, .5, .75, .975])] for x in scenario_results]

AR_convergence_check = [mean(x[:n_human_recovered_reps][1:i,n_steps] .+ x[:n_human_infections_reps][1:i,n_steps]) for x in scenario_results, i in 1:n_reps]

plot(transpose(AR_convergence_check), label=label=reshape(these_labs, 1, length(scenario_results)), xlabel="Repititions", ylabel="Mean Attack Rate", palette = :Dark2_8, linewidth=2)
plot!(size=(800,600))
png(string("plots/AR_conv", savename(this_sim_dict), ".png"))



mean_bites_per_person = [mean(x[4]) for x in scenario_results]

max_bites_per_person = [maximum(x[4]) for x in scenario_results]


