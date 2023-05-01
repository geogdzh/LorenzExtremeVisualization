# make holding time plot, generator, and network graph for static ρ=28 scenario 
using MarkovChainHammer: generator
include("./analyze_util.jl")

# here I need one static run, and then two static and one delta
markov_chain, dt = read_markov_chain("./data/lorenz28.hdf5")
mc_ag = aggregate(markov_chain)

## generator shown in paper:
Q_28 = generator(markov_chain; dt=dt)

Q_condensed = generator(mc_ag; dt=dt)

using Latexify
latexify(round.(Q_condensed, digits=3))

## network graph

include("../visualize/graph_plot.jl")