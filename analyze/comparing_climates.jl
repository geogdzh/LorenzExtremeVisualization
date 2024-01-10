# generate butterflies and holding time comparison
using HDF5
using GLMakie
using MarkovChainHammer.TransitionMatrix: holding_times, generator
include("./analyze_util.jl")

##

mc28, dt, x = read_markov_chain("./data/lorenz28.hdf5", include_x=true)
mc28_ag = aggregate(mc28)
z = x[3,:]
nsteps = 2100
tlist = collect(0:dt:dt*(length(x)-1))
color_choices = ["blue", "violet", "red"]
colorlist = [color_choices[mc28_ag[i]] for i in 1:nsteps]

fig = Figure(resolution=(1000, 500))
ax = Axis(fig[1, 1]; xlabel="Time", ylabel="z")
scatter!(ax, tlist[1:nsteps], z[1:nsteps], color=colorlist)
elem_1 = [MarkerElement(color=:blue, marker=:circle, markersize=12)]
elem_2 = [MarkerElement(color=:violet, marker=:circle, markersize=12)]
elem_3 = [MarkerElement(color=:red, marker=:circle, markersize=12)]
axislegend(ax, [elem_1, elem_2, elem_3], ["A", "B", "C"], position=:lt)
save("figs/z_trajectory.png", fig)
fig

##

include("../visualize/partitioned_butterfly.jl")

##

include("../visualize/graph_plot.jl")
