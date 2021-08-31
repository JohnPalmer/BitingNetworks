these_labs = [string(x) for x in keys(human_distribution)]


maximum(scenario_results[1][1])

using Plots.PlotMeasures
these_plots = []

this_max = maximum([maximum(x[1]) for x in scenario_results])/n_humans

steps_to_plot = 15

for i in 1:length(human_distributions)
  n_human_infections_reps = scenario_results[i][1]./n_humans
  n_human_infections_reps_median = transpose(median(n_human_infections_reps, dims=1))

  this_p = plot(1:steps_to_plot, transpose(n_human_infections_reps[ :, 1:steps_to_plot]), w=.2, color=:grey, legend=:none, yaxis = ("Infected", (0, this_max), 0:.5:this_max), xaxis = ("Time"), margin = 1cm)


  plot!(this_p, 1:steps_to_plot,n_human_infections_reps_median[1:steps_to_plot], w=2, color=:black, yaxis = ("Infected", (0, 1), 0:.2:1))

  push!(these_plots, this_p)
end

max_bites = maximum([maximum(x[4]) for x in scenario_results])

kd_max = maximum(kde(vec(scenario_results[1][4])).density)

these_bite_plots = [density(scenario_results[i][:human_bite_distribution], xaxis = ("N Bites", (0, max_bites), 0:ceil(max_bites/10):max_bites), yaxis = ("Density", (-0.01, kd_max)), legend=:none, normed=true, linewidth=2, title=these_labs[i], margin = 1cm) for i in 1:length(scenario_results)]

for i in 1:length(these_plots)
  this_p = plot( 
    these_bite_plots[i], these_plots[i], layout=(1,2))
  plot!(this_p, size=(600,400))
  savefig(this_p, string("plots/epi_comp_", savename(this_sim_dict_a), human_distribution_names[i], ".html"))
end
