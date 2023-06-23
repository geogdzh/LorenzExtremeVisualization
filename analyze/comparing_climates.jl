# generate butterflies and holding time comparison

include("./analyze_util.jl")

# here I need two static runs and one delta one
mc26, dt = read_markov_chain("./data/lorenz26.hdf5")
mc32, dt = read_markov_chain("./data/lorenz32.hdf5")
mcdelta, dt_delta = read_markov_chain("./data/lorenz-changing-10e6.hdf5")

##

include("../visualize/partitioned_butterfly.jl")

##
# now, holding time comparison: (aggregated)

mc26_ag = aggregate(mc26)
mc32_ag = aggregate(mc32)
mcdelta_ag = aggregate(mcdelta)

ht_delta = holding_times(mcdelta_ag, 3; dt=dt)
ht_26 = holding_times(mc26_ag, 3; dt=dt)
ht_32 = holding_times(mc32_ag, 3; dt=dt)

begin
    fig = Figure(resolution=(1600, 1200))
    lims = [(0,1.2 ), (0,2), (0,0.2)]
    for i in 1:3
        start, stop = lims[i]
        bin_num=15
        bin_width = (stop-start)/bin_num
        ax = Axis(fig[i,1], title="$i (ρ=26)") 
        hist!(ax, ht_26[i], normalization=:probability,bins=start:bin_width:stop)
        xlims!(start, stop)
        # ylims!(0, 0.4)
        #
        ax = Axis(fig[i,2], title="$i (ρ changing)") 
        hist!(ax, ht_delta[i], normalization=:probability,bins=bins=start:bin_width:stop)
        xlims!(start, stop)
        # ylims!(0, 0.4)
        #
        ax = Axis(fig[i,3], title="$i (ρ=32)") 
        hist!(ax, ht_32[i], normalization=:probability,bins=bins=start:bin_width:stop)
        xlims!(start, stop)
        # ylims!(0, 0.4)
    end
    fig
    save("figs/holding_time_histogram.png", fig)
end
