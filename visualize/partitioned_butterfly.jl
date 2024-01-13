using HDF5, GLMakie, LinearAlgebra, Statistics

include("../analyze/analyze_util.jl")

mc26, dt, x26 = read_markov_chain("./data/lorenz26.hdf5"; include_x=true)
mc32, dt, x32 = read_markov_chain("./data/lorenz32.hdf5"; include_x=true)
# mcdelta, dt_delta, xdelta = read_markov_chain("./data/lorenz-changing-10e7.hdf5"; include_x=true) 


# create colors for the plot
# non-custom, see https://docs.juliaplots.org/latest/generated/colorschemes/
# color_choices = cgrad(:rainbow, 12, categorical=true) # 12 is the number of states
# colormap = :Paired_12
color_choices = [:yellow, :darkgrey, :lightgrey, :magenta ,:blue, :cyan, :lightgreen, :green ,:lightblue, :red , :purple, :orange]
# color_choices = [cgrad(colormap)[i] for i in 1:12]

colors26 = []
for i in 1000:Int(10e6) #eachindex(mc26)
    push!(colors26, color_choices[mc26[i]])
end
# tuple_timeseries_26 = [(x26[1, i], x26[2, i], x26[3, i]) for i in 1:size(x26, 2)]
tuple_timeseries_26 = [(x26[1, i], x26[2, i], x26[3, i]) for i in 1000:Int(10e6)]

colors32 = []
for i in 1000:Int(10e6)#eachindex(mc32)
    push!(colors32, color_choices[mc32[i]])
end
# tuple_timeseries_32 = [(x32[1, i], x32[2, i]+60, x32[3, i]) for i in 1:size(x32, 2)]
tuple_timeseries_32 = [(x32[1, i], x32[2, i]+60, x32[3, i]) for i in 1000:Int(10e6)]

# colorsdelta = []
# for i in eachindex(mcdelta)
#     push!(colorsdelta, color_choices[mcdelta[i]])
# end
# tuple_timeseries_delta = [(xdelta[1, i], xdelta[2, i]+40, xdelta[3, i]) for i in 1:size(xdelta, 2)]


# plotting
fig = Figure(resolution=(2100, 1000))
ax = LScene(fig[1:2, 1:2]; show_axis=true) # show axis!
lines!(ax, tuple_timeseries_26, color=colors26)
lines!(ax, tuple_timeseries_32, color=colors32)
# lines!(ax, tuple_timeseries_delta, color=colorsdelta)
rotate_cam!(ax.scene, (pi/6, -π / 4, 0))
display(fig)

# to record video:
# last_time_index = minimum([60 * 15 * 2, length(timeseries)])
# time_indices = 1:last_time_index

# function change_function(time_index)
#     phase = 2π / (60 * 15)
#     rotate_cam!(ax.scene, (0, phase, 0))
# end

#record(change_function, fig, "lorenz_animation_attractor_2.mp4", time_indices; framerate=60)


