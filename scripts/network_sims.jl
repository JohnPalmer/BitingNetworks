using Random, Distributions
using ProgressMeter, Plots, StatsPlots
gr()

Threads.nthreads()

Random.seed!(123)

n_steps = 100
n_reps = 300
n_humans = 100
n_mosquitoes = 400

neg_bin_r = 2
neg_bin_p = .4

transmission_prob = .05

human_distribution = NegativeBinomial(neg_bin_r, neg_bin_p)

# human_distribution = Beta(2, 5)
# mosquito_distribution = Beta(.2, 2)

# human_distribution = Uniform(0, 1)
# mosquito_distribution = Uniform(0, 1)

# human_probs = rand(human_distribution, n_humans)
# mosquito_probs = rand(mosquito_distribution, n_mosquitoes)

mosquito_probs1 = fill(1/n_humans, (n_mosquitoes,1))

human_probs1 = rand(human_distribution, n_humans)/(sum(mosquito_probs))

human_probs2 = fill(mean(human_probs1), (n_humans,1))
mosquito_probs2 = fill(mean(mosquito_probs1), (n_mosquitoes,1))

scenario_results = []

scenarios = [[mosquito_probs1, human_probs1],[mosquito_probs2, human_probs2] ]

scenario = scenarios[2]
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

  p = Progress(n_reps)

  Threads.@threads for r = 1:n_reps

  # creating first infection
  status_humans = zeros(Int8, n_humans)
  status_humans[rand(1:n_humans, 1)[1]] = 1

  status_mosquitoes = zeros(Int8, n_mosquitoes)

  n_human_infections = Vector{Int64}(undef, n_steps)
  n_mosquito_infections = Vector{Int64}(undef, n_steps)

  n_human_infections[1] = 1
  n_mosquito_infections[1] = 0

    for s = 2:n_steps

      i_hs = findall(status_humans .== 1)
      s_ms = findall(status_mosquitoes .==0)
      these_iterations = 1:length(i_hs)*length(s_ms)

    # human-to-mosquito infections
      for i in these_iterations
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
      for i in these_iterations
        h = s_hs[cld(i, length(i_ms))]
        m = i_ms[mod1(i, length(i_ms))]
        if rand(Bernoulli(human_probs[h]*mosquito_probs[m]*transmission_prob), 1)[1]
          status_humans[h] = 1
        end
      end

      n_human_infections[s] = sum(status_humans)

    end

    n_human_infections_reps[r, :] = n_human_infections

    n_mosquito_infections_reps[r, :] = n_mosquito_infections

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

scenario_result = [n_human_infections_reps, human_R0, mosquito_R0, human_bite_distribution, mosquito_bite_distribution]

push!(scenario_results, scenario_result)

end

scenario_results[1][2]
scenario_results[2][2]

mean(scenario_results[1][4])
mean(scenario_results[2][4])

these_plots = []

for i in 1:2
  n_human_infections_reps = scenario_results[i][1]
  n_human_infections_reps_median = transpose(median(n_human_infections_reps, dims=1))

  this_p = plot(1:n_steps, transpose(n_human_infections_reps), w=.2, color=:grey, legend=:none)

  plot!(this_p, 1:n_steps,n_human_infections_reps_median, w=2, color=:black)

  push!(these_plots, this_p)
end

plot(these_plots[1], these_plots[2])