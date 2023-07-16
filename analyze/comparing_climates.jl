# generate butterflies and holding time comparison
using HDF5
using GLMakie
using MarkovChainHammer.TransitionMatrix: holding_times
include("./analyze_util.jl")

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
    xlims = [(0,1.2 ), (0,1.99), (0,0.19)]
    ylims = [(0,0.3), (0,0.7), (0,0.4)]
    for i in 1:3
        start, stop = xlims[i]
        bin_num=15
        bin_width = (stop-start)/bin_num
        kwargs = (; xlabel=(i==3 ? "Holding time" : ""), titlesize=22, ylabelsize=20, xlabelsize=20, yticklabelsize=18, xticklabelsize=18)
        ax = Axis(fig[i,1], title="State $i (ρ=26)", ylabel="Probability"; kwargs...) 
        hist!(ax, ht_26[i], normalization=:probability,bins=start:bin_width:stop)
        xlims!(start, stop)
        ylims!(ylims[i]...)
        #
        ax = Axis(fig[i,2], title="State $i (ρ changing)"; kwargs...) 
        hist!(ax, ht_delta[i], normalization=:probability,bins=bins=start:bin_width:stop)
        xlims!(start, stop)
        ylims!(ylims[i]...)
        hideydecorations!(ax, grid=false)
        #
        ax = Axis(fig[i,3], title="State $i (ρ=32)"; kwargs...) 
        hist!(ax, ht_32[i], normalization=:probability,bins=bins=start:bin_width:stop)
        xlims!(start, stop)
        ylims!(ylims[i]...)
        hideydecorations!(ax, grid=false)
    end
    fig
    save("figs/newest/holding_time_histogram.png", fig)
end
