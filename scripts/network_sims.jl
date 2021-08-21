using Random, Distributions
using ProgressMeter, Plots, StatsPlots
gr()

Threads.nthreads()

Random.seed!(123)

n_steps = 50
n_reps = 300
n_humans = 1000
n_mosquitoes = 4000

human_bite_mean = 10

transmission_prob = .05

human_distribution = Poisson(human_bite_mean)

# human_distribution = Beta(2, 5)
# mosquito_distribution = Beta(.2, 2)

# human_distribution = Uniform(0, 1)
# mosquito_distribution = Uniform(0, 1)

# human_probs = rand(human_distribution, n_humans)
# mosquito_probs = rand(mosquito_distribution, n_mosquitoes)

# human_probs = fill(.5, (n_humans,1))
# mosquito_probs = fill(.5, (n_mosquitoes,1))



mosquito_probs = fill(1/(n_humans*human_bite_mean), (n_mosquitoes,1))

human_probs = rand(human_distribution, n_humans)/sum(mosquito_probs)

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

density(human_bite_distribution)

density(mosquito_bite_distribution)

n_human_infections_reps = Array{Int64}(undef, n_reps, n_steps)

n_mosquito_infections_reps = Array{Int64}(undef, n_reps, n_steps)

@showprogress Threads.@threads for r = 1:n_reps

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
end

n_human_infections_reps_median = transpose(median(n_human_infections_reps, dims=1))

plot(1:n_steps, transpose(n_human_infections_reps), w=.2, color=:grey, legend=:none)

plot!(1:n_steps,n_human_infections_reps_median, w=2, color=:black)

R0s = Vector(undef, n_reps)

for i in 1:n_reps
  x = n_human_infections_reps[i,:]
  if maximum(x) > 1
    R0s[i] = minimum(x[x.>1])
  else
    R0s[i] = 0
  end
end

mean(R0s)