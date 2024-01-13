using HDF5
using GLMakie
using MarkovChainHammer.TransitionMatrix: holding_times, generator
include("./analyze_util.jl")

##
# using basic holding time runs:

# mc26, dt = read_markov_chain("./data/lorenz26.hdf5")
# mc32, dt = read_markov_chain("./data/lorenz32.hdf5")
# mcdelta, dt_delta = read_markov_chain("./data/lorenz-changing-10e6.hdf5")

# mc26_ag = aggregate(mc26)
# mc32_ag = aggregate(mc32)
# mcdelta_ag = aggregate(mcdelta)

# ht_delta = holding_times(mcdelta_ag, 3; dt=dt)
# ht_26 = holding_times(mc26_ag, 3; dt=dt)
# ht_32 = holding_times(mc32_ag, 3; dt=dt)


##
# using instead the extra holding times:
hfile = h5open(pwd() * "/data/lorenz26_holding-times.hdf5")
ht_26 = [read(hfile["ht1"]), read(hfile["ht2"]), read(hfile["ht3"])]
close(hfile)

hfile = h5open(pwd() * "/data/lorenz32_holding-times.hdf5")
ht_32 = [read(hfile["ht1"]), read(hfile["ht2"]), read(hfile["ht3"])]
close(hfile)

hfile = h5open(pwd() * "/data/ensemble-holding-times.hdf5")
ht_delta = [read(hfile["ht1"]), read(hfile["ht2"]), read(hfile["ht3"])]
close(hfile)

begin
    fig = Figure(resolution=(1600, 1200))
    xlims = [(0,1.2), (0,0.6), (0,0.2)]
    ylims = [(0,0.15), (0,0.3), (0,0.15)]
    states = ["A", "B", "C"]
    xlabels = [([0.0, 0.5, 1.0], ["0.0", "0.5", "1.0"]), ([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5"]), ([0.0, 0.05, 0.1, 0.15], ["0.0", "0.05", "0.1", "0.15"])]
    for i in 1:3
        start, stop = xlims[i]
        bin_num= 39 #(i==2 ? 40 : 40)
        bin_width = (stop-start)/bin_num
        kwargs = (; xlabel=(i==3 ? "Holding time" : ""), titlesize=22, ylabelsize=20, xlabelsize=20, yticklabelsize=18, xticklabelsize=18, xticks = xlabels[i])
        ax = Axis(fig[i,1], title="Macrostate $(states[i]) (ρ=26)", ylabel="Probability"; kwargs...) 
        hist!(ax, ht_26[i], normalization=:probability,bins=start:bin_width:stop)
        xlims!(start, stop)
        ylims!(ylims[i]...)
        #
        ax = Axis(fig[i,2], title="Macrostate $(states[i]) (ρ changing)"; kwargs...) 
        hist!(ax, ht_delta[i], normalization=:probability,bins=bins=start:bin_width:stop)
        xlims!(start, stop)
        ylims!(ylims[i]...)
        hideydecorations!(ax, grid=false)
        #
        ax = Axis(fig[i,3], title="Macrostate $(states[i]) (ρ=32)"; kwargs...) 
        hist!(ax, ht_32[i], normalization=:probability,bins=bins=start:bin_width:stop)
        xlims!(start, stop)
        ylims!(ylims[i]...)
        hideydecorations!(ax, grid=false)
    end
    fig
    save("figs/holding_time_histogram.png", fig)
end