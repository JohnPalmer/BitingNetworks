using DrWatson
quickactivate(@__DIR__, "BitingNetworks")
using Bites, Random, Distributions
using ProgressMeter, Plots, StatsPlots, KernelDensity, DataFrames, CSV
plotlyjs()

Threads.nthreads()

Random.seed!(123)

n_steps = 50
n_reps = 1000
n_humans = 1000
n_mosquitoes = 4000

transmission_prob = .15

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
  mixed_norms = vcat(subpop_a, subpop_b),
  tlevy1p1 = Truncated(Levy(1.1, .0001), 1, 1000),
  tlevy2p1 = Truncated(Levy(2.1, .0001), 1, 1000),
  tlevy3p1 = Truncated(Levy(3.1, .0001), 1, 1000),
  tlevy4p1 = Truncated(Levy(4.1, .0001), 1, 1000),
  tlevy5p1 = Truncated(Levy(5.1, .0001), 1, 1000),
  tlevy6p1 = Truncated(Levy(6.1, .0001), 1, 1000),
)

this_set_name = "cueml16"

human_distribution_names = join([string(x) for x in keys(human_distributions)], ".")

this_sim_dict_a = @strdict n_steps n_reps n_humans n_mosquitoes transmission_prob expected_bites human_infection_time mosquito_life_span

this_sim_dict = @strdict n_steps n_reps n_humans n_mosquitoes transmission_prob expected_bites human_infection_time mosquito_life_span this_set_name

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

@tagsave(datadir("sims", savename(scenario_result, "jld2")), tostringdict(ntuple2dict(scenario_result)), safe=true)

push!(scenario_results, scenario_result)

end

## SUMMARY ANALYSIS ##

# main labels
these_labs = [string(x) for x in keys(human_distributions)]

# Bootstrap R0
R0_boot_reps = [mean(rand(x.human_R0_reps, length(x.human_R0_reps))) for i in 1:1000, x in scenario_results]

R0_boot = DataFrame(R0_boot_reps, these_labs) 

R0_boot_long = stack(R0_boot, 1:length(these_labs))
this_p = @df R0_boot_long boxplot(:variable, :value, ylabel="R0", legend=:none)
png(this_p, string("plots/R0_boot_", savename(this_sim_dict), ".png"))

# R0 Convergence
human_R0_convergence_checks = [x[:human_R0_converge_check] for x in scenario_results]

Plots.plot(human_R0_convergence_checks, label=reshape(these_labs, 1, length(scenario_results)), xlabel="Repetitions", ylabel="Mean R0", legend=:bottomright, linewidth=2, palette = :Dark2_8)
Plots.plot!(size=(800,600))
png(string("plots/R0_conv", savename(this_sim_dict), ".png"))

# AR Convergence
AR_convergence_check = [mean(x[:n_human_recovered_reps][1:i,n_steps] .+ x[:n_human_infections_reps][1:i,n_steps])/n_humans for x in scenario_results, i in 1:n_reps]

Plots.plot(transpose(AR_convergence_check), label=label=reshape(these_labs, 1, length(scenario_results)), xlabel="Repititions", ylabel="Mean Attack Rate", palette = :Dark2_8, linewidth=2)
Plots.plot!(size=(800,600))
png(string("plots/AR_conv", savename(this_sim_dict), ".png"))

# Point estimate summary

summary_estimates = DataFrame(

human_distribution = these_labs,

R0 = [x.human_R0 for x in scenario_results],

R0_ul90 = [quantile(R0_boot_reps[:, i], .95) for i in 1:length(scenario_results)],

R0_ll90 = [quantile(R0_boot_reps[:, i], .05) for i in 1:length(scenario_results)],

AR = [mean(x[:n_human_recovered_reps][:,n_steps] .+ x[:n_human_infections_reps][:,n_steps])/n_humans for x in scenario_results],

AR_ul90 = [quantile(x[:n_human_recovered_reps][:,n_steps] .+ x[:n_human_infections_reps][:,n_steps], .95)/n_humans for x in scenario_results],

AR_ll90 = [quantile(x[:n_human_recovered_reps][:,n_steps] .+ x[:n_human_infections_reps][:,n_steps], .05)/n_humans for x in scenario_results],


bites_per_person_mean = [mean(x.human_bite_distribution) for x in scenario_results],

bites_per_person_max = [maximum(x.human_bite_distribution) for x in scenario_results],

bites_per_person_min = [minimum(x.human_bite_distribution) for x in scenario_results],

bites_per_person_ul90 = [quantile(vec(x.human_bite_distribution), .95) for x in scenario_results],

bites_per_person_ll90 = [quantile(vec(x.human_bite_distribution), .05) for x in scenario_results],

bites_per_person_ul50 = [quantile(vec(x.human_bite_distribution), .75) for x in scenario_results],

bites_per_person_ll50 = [quantile(vec(x.human_bite_distribution), .25) for x in scenario_results],

bites_per_mosquito_mean = [mean(x.mosquito_bite_distribution) for x in scenario_results],

bites_per_mosquito_max = [maximum(x.mosquito_bite_distribution) for x in scenario_results],

bites_per_mosquito_min = [minimum(x.mosquito_bite_distribution) for x in scenario_results],

bites_per_mosquito_ul90 = [quantile(vec(x.mosquito_bite_distribution), .95) for x in scenario_results],

bites_per_mosquito_ll90 = [quantile(vec(x.mosquito_bite_distribution), .05) for x in scenario_results],

bites_per_mosquito_ul50 = [quantile(vec(x.mosquito_bite_distribution), .75) for x in scenario_results],

bites_per_mosquito_ll50 = [quantile(vec(x.mosquito_bite_distribution), .25) for x in scenario_results]
)


CSV.write(datadir("sim_summaries", savename(this_sim_dict, "csv")), summary_estimates)

@tagsave(datadir("sim_summaries", savename(this_sim_dict, "jld2")), tostringdict(struct2dict(summary_estimates)), safe=true)