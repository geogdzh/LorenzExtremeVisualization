# generate butterflies and holding time comparison

include("./lorenz_embedding.jl")

# here I need two static runs and one delta one
mc26, dt = get_markov_chain("./data/lorenz26.hdf5")
mc32, dt = get_markov_chain("./data/lorenz32.hdf5")
mcdelta, dt_delta = get_markov_chain("./data/lorenz-changing-10e6.hdf5") # is this the one we want to use??
##

function partitioned_butterfly(filename)
    @info "opening data"
    hfile = h5open(filename)
    x = read(hfile["x"])
    dt = read(hfile["dt"]) #does this even work
    ρ = read(hfile["rho"])
    close(hfile)
    @info "applying embedding"
    markov_indices = zeros(Int, size(x,2))
    for i in ProgressBar(eachindex(markov_indices))
        markov_indices[i] = embedding(x[:, i])
    end

    # create colors for the plot
    colors = []
    # non-custom, see https://docs.juliaplots.org/latest/generated/colorschemes/
    # color_choices = cgrad(:rainbow, 12, categorical=true) # 12 is the number of states
    color_choices = [:yellow, :darkgrey, :lightgrey, :pink ,:blue, :lightblue, :lightgreen, :green ,:red, :magenta , :purple, :orange]

    for i in eachindex(x)
        push!(colors, color_choices[markov_indices[i]]) #no clue why this is throwing an error
    end

    tuple_timeseries = Tuple.(x)

    # everything is done for plotting
    fig = Figure(resolution=(1000, 700))
    ax = LScene(fig[1:2, 1:2]; show_axis=false)
    lines!(ax, tuple_timeseries, color=colors)
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

end

# generate figures
partitioned_butterfly("./data/lorenz26.hdf5") # need to make sure the data locations are consistent!!
partitioned_butterfly("./data/lorenz32.hdf5")
partitioned_butterfly("./data/lorenz-changing-10e6.hdf5")

##
# now, holding time comparison: (aggregated)

mc26-ag, dt = get_markov_chain("./data/lorenz26.hdf5"; ag=true)
mc32-ag, dt = get_markov_chain("./data/lorenz32.hdf5")
mcdelta-ag, dt = get_markov_chain("./data/lorenz-changing-10e6.hdf5") # is this the one we want to use??

ht_delta = holding_times(mc-delta-ag, 3; dt=dt)
ht_26 = holding_times(mc26-ag, 3; dt=dt)
ht_32 = holding_times(mc32-ag, 3; dt=dt)

begin
    fig = Figure(resolution=(1600, 1200))
    for i in 1:3
        ax = Axis(fig[i,1], title="$i (ρ=26)") 
        hist!(ax, ht_26[i], normalization=:probability,nbins=10)
    end
    for i in 1:3
        ax = Axis(fig[i,2], title="$i (ρ changing)") 
        hist!(ax, ht_delta[i], normalization=:probability,nbins=10)
    end
    for i in 1:3
        ax = Axis(fig[i,3], title="$i (ρ=32)") 
        hist!(ax, ht_32[i], normalization=:probability,nbins=10)
    end
    # fig[0, 1] = Label(fig, "ρ=26")
    # xlim!(fig[2,3], (0,4))
    fig
end


