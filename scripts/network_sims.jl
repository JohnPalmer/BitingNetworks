using Random, Distributions
using ProgressMeter, Plots, StatsPlots
gr()

Threads.nthreads()

Random.seed!(123)

n_steps = 50
n_reps = 15
n_humans = 1000
n_mosquitoes = 10000

# how long before infected human recovers
human_infection_time = 3

# how long do mosquitoes live. Since mosquitoes do not usually get accute arbovirus infections, the life span is the key: They will go from susceptible->infected->dead/reborn. (With assumption of population stability, deaths simply put the agent back to susceptible status.) "Unlike arboviral infections in humans, which are usually acute, arboviral infections in mosquitoes are persistent. Once the infection is established, the mosquito remains infected for the rest of its life." https://www.sciencedirect.com/science/article/pii/S1931312819303701
mosquito_life_span = 20

neg_bin_r = 2
neg_bin_p = .1

neg_bin_r_a = 5
neg_bin_p_a = .2

neg_bin_r_b = 10
neg_bin_p_b = .35

transmission_prob = .1

human_distribution = NegativeBinomial(neg_bin_r, neg_bin_p)

human_distribution_a = NegativeBinomial(neg_bin_r_a, neg_bin_p_a)

human_distribution_b = NegativeBinomial(neg_bin_r_b, neg_bin_p_b)

# human_distribution = Beta(2, 5)
# mosquito_distribution = Beta(.2, 2)
mosquito_probs1 = fill(1/n_humans, (n_mosquitoes,1))

human_probs1 = rand(human_distribution, n_humans)/(sum(mosquito_probs1))

human_probs1a = rand(human_distribution_a, n_humans)/(sum(mosquito_probs1))

human_probs1b = rand(human_distribution_b, n_humans)/(sum(mosquito_probs1))

human_probs2 = fill(mean(human_probs1), (n_humans,1))
mosquito_probs2 = fill(mean(mosquito_probs1), (n_mosquitoes,1))


human_distribution3 = Uniform(0, mean(human_probs1)*2)
mosquito_distribution3 = Uniform(0, mean(mosquito_probs1)*2)

human_probs3 = rand(human_distribution3, n_humans)
mosquito_probs3 = rand(mosquito_distribution3, n_mosquitoes)


scenario_results = []

scenarios = [[mosquito_probs1, human_probs1],[mosquito_probs1, human_probs1a], [mosquito_probs1, human_probs1b],[mosquito_probs2, human_probs2], [mosquito_probs3, human_probs3] ]

scenario = scenarios[3]

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

  p = Progress(n_reps)

  Threads.@threads for r = 1:n_reps

  # vectors of infections statusas follows: 0 = susceptible, >0 = infected, <0 = recovered. Everyone starts susceptible
  status_humans = zeros(Int8, n_humans)
  status_mosquitoes = zeros(Int8, n_mosquitoes)

  # vector of mosquito ages
  age_mosquitoes = rand(DiscreteUniform(0, mosquito_life_span), n_mosquitoes)

  # creating first infection
  status_humans[rand(1:n_humans, 1)[1]] = 1

  n_human_infections = Vector{Int64}(undef, n_steps)
  n_mosquito_infections = Vector{Int64}(undef, n_steps)

  n_human_recovered = Vector{Int64}(undef, n_steps)

  n_human_infections[1] = 1
  n_mosquito_infections[1] = 0
  n_human_recovered[1] = 0

    for s = 2:n_steps

      # update statuses based on time n_steps
      status_humans[findall(status_humans .> 0)] .+= 1
      status_humans[findall(status_humans .> human_infection_time)] .= -1

      # all mosquitoes age by 1 step
      age_mosquitoes .+= 1
      # mosquitoes over max age die and are reborn (under assumption of stable population) as susceptible
      status_mosquitoes[findall(age_mosquitoes .> mosquito_life_span)] .= 0
      age_mosquitoes[findall(age_mosquitoes .> mosquito_life_span)] .= 0

      # find indexes of infected humans (status>0)
      i_hs = findall(status_humans .> 0)
      # find indexes of susecptible mosquitoes (status == 0)
      s_ms = findall(status_mosquitoes .==0)
      
      # calculate number of iterations needed to cycle through all combinations. (This is useful if parallizing over these iterations, since it allows for more efficient distribution of the tasks)
      these_iterations = 1:length(i_hs)*length(s_ms)

    # human-to-mosquito infections
      for i in these_iterations
        h = i_hs[cld(i, length(s_ms))]
        m = s_ms[mod1(i, length(s_ms))]
        if rand(Bernoulli(human_probs[h]*mosquito_probs[m]*transmission_prob), 1)[1]
          status_mosquitoes[m] = 1
        end
      end

    n_mosquito_infections[s] = sum(status_mosquitoes .> 0)

    s_hs = findall(status_humans .== 0)
    i_ms = findall(status_mosquitoes .==1)
    these_iterations = 1:length(s_hs)*length(i_ms)

      # mosquito-human infections
      for i in these_iterations
        h = s_hs[cld(i, length(i_ms))]
        m = i_ms[mod1(i, length(i_ms))]
        if rand(Bernoulli(human_probs[h]*mosquito_probs[m]*transmission_prob), 1)[1]
          status_humans[h] = 1
        end
      end

      n_human_infections[s] = sum(status_humans.>0)
      n_human_recovered[s] = sum(status_humans.<0)

    end

    n_human_infections_reps[r, :] = n_human_infections

    n_mosquito_infections_reps[r, :] = n_mosquito_infections

    n_human_recovered_reps[r, :] = n_human_recovered

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

scenario_result = [n_human_infections_reps, human_R0, mosquito_R0, human_bite_distribution, mosquito_bite_distribution, n_human_recovered_reps]

push!(scenario_results, scenario_result)

end

R0 = [x[2] for x in scenario_results]

scenario_results[1][6][:,n_steps] .+ scenario_results[1][1][:,n_steps] 

AR = [ [mean(x[6][:,n_steps] .+ x[1][:,n_steps]), quantile(x[6][:,n_steps] .+ x[1][:,n_steps], [.025, .25, .5, .75, .975])]  for x in scenario_results]

mean_bites_per_person = [mean(x[4]) for x in scenario_results]


these_plots = []

this_max = maximum([maximum(x[1]) for x in scenario_results])

for i in 1:length(scenarios)
  n_human_infections_reps = scenario_results[i][1]
  n_human_infections_reps_median = transpose(median(n_human_infections_reps, dims=1))

  this_p = plot(1:n_steps, transpose(n_human_infections_reps), w=.2, color=:grey, legend=:none, yaxis = ("Infected", (0, this_max)), xaxis = ("Time"))

  plot!(this_p, 1:n_steps,n_human_infections_reps_median, w=2, color=:black)

  push!(these_plots, this_p)
end

max_bites = maximum([maximum(x[4]) for x in scenario_results])

these_bite_plots = [density(x[4], xaxis = ("N Bites", (0, max_bites), 0:40:max_bites), yaxis = ("Density", (0, .1), 0:.02:.1), legend=:none) for x in scenario_results]

l = @layout [a b c d e; f g h i j]

plot(these_plots[1], these_plots[2], these_plots[3],these_plots[4], these_plots[5], these_bite_plots[1], these_bite_plots[2], these_bite_plots[3], these_bite_plots[4], these_bite_plots[5], layout = l)