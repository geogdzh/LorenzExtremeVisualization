using HDF5, GLMakie, LinearAlgebra, Statistics, Random
# using MarkovChainHammer.BayesianMatrix
# using MarkovChainHammer.TransitionMatrix: steady_state
# using MarkovChainHammer.TransitionMatrix: perron_frobenius
using MarkovChainHammer
using SparseArrays, Graphs, GraphMakie, Printf
##
include("../analyze/analyze_util.jl")

mc, dt, x = read_markov_chain("data/lorenz26.hdf5"; include_x=true)
mc26 = aggregate(mc[1:10000000])
mc, dt, x = read_markov_chain("data/lorenz32.hdf5"; include_x=true)
mc32 = aggregate(mc[1:10000000])


gen26 = generator(mc26; dt=dt)
gen32 = generator(mc32; dt=dt)

# generate numerical matrices being shown in paper
using Latexify, LaTeXStrings
latexify(round.(gen26, digits=2))
e,r = exit_rate(gen26)
latexify(round.(e, digits=2))
latexify(round.(1 ./ r, digits=2))
latexify(round.(gen32, digits=2))
e,r = exit_rate(gen32)
latexify(round.(e, digits=2))
latexify(round.(1 ./ r, digits=2))

## generate two side by side graph network plots, for each of the matrices

function edge_graph(ax, mc)
    Q = mean(BayesianGenerator(mc; dt = dt))
    # Generator
    g_Q = DiGraph(Q')
    Q_prim = zeros(size(Q))
    for i in 1:size(Q)[1]
        Q_prim[:, i] .= -Q[:, i] / Q[i, i]
        Q_prim[i, i] = 0.#-1 / Q[i, i] #cahnged to zero to hide self-loops
    end
    hts = [round(-1 / Q[i, i], digits=3) for i in 1:3]

    elabels = string.([round(Q_prim[i]; digits=2) for i in 1:ne(g_Q)])
    # [@sprintf("%.0e", node_labels[i]) for i in 1:nv(G)]
    # transparancy = [Q_prim[i] for i in 1:ne(g_Q)]
    transparancy = [(Q_prim[i] > eps(10^6 * 1.0) ? maximum([Q_prim[i], 0.2]) : 0.0)  for i in 1:ne(g_Q)]
    elabels_color = [(:black, transparancy[i] > eps(100.0)) for i in 1:ne(g_Q)]
    # edge_color_Q = [(cgrad(colormap)[(i - 1)÷12 + 1], transparancy[i]) for i in 1:ne(g_Q)]
    edge_color_Q = [(colors[(i - 1)÷3 + 1], transparancy[i]) for i in 1:ne(g_Q)]
    # node_color = [(cgrad(colormap)[i]) for i in 1:nv(g_Q)]
    node_color = [(colors[i]) for i in 1:nv(g_Q)]
    edge_attr = (; linestyle=[:dot, :dash, :dash, :dash, :dot, :dash, :dash, :dash, :dot])
    elabels_fontsize = 36
    nlabels_fontsize = 36
    node_size = 50.0

    edge_width_Q = [10.0 for i in 1:ne(g_Q)]
    arrow_size_Q = [40.0 for i in 1:ne(g_Q)]
    # node_labels_Q = [L"A (normal)   \newline $<T_A>=%$(hts[1])$   ", L"B (medium)   \n $<T_B>=%$(hts[2])$   ", L"C (high)  \n $<T_C>=%$(hts[3])$  "]#repr.(1:nv(g_Q))
    node_labels_Q = [L"<T_A>=%$(hts[1])         ", L"$<T_B>=%$(hts[2])$   ", L"$<T_C>=%$(hts[3])$  "]#repr.(1:nv(g_Q))

    nlabels_align = [(:right,:bottom), (:right,:center), (:right,:bottom)]

    #                   ; edge_color=edge_color_Q, edge_width=edge_width_Q
    kwargs_edges = (; elabels=elabels, elabels_color=elabels_color, elabels_fontsize=elabels_fontsize, edge_color=edge_color_Q, edge_width=edge_width_Q)
    kwargs_nodes  = (; node_color=node_color, node_size=node_size, nlabels=node_labels_Q, nlabels_fontsize=nlabels_fontsize, nlabels_align=nlabels_align)
    kwargs_arrows = (; arrow_size=arrow_size_Q)
    # kwargs_edges = (; elabels_color=elabels_color, edge_color=edge_color_Q, edge_width=edge_width_Q)
    # kwargs_nodes  = (; node_color=node_color, node_size=node_size)
    # kwargs_arrows = (; arrow_size=arrow_size_Q)
        
    # p_Q = graphplot!(ax, g_Q; kwargs_edges..., kwargs_nodes..., kwargs_arrows...)
    graphplot!(ax, g_Q; kwargs_edges..., kwargs_nodes..., kwargs_arrows..., elabels_distance=20)#, layout = Stress())
    hidedecorations!(ax)
    hidespines!(ax)
end
##

fig = Figure(resolution=(1500, 560))
colors = ["blue", "violet", "red"]

ax_26 = Axis(fig[1, 1]; title="ρ=26", titlesize=36)
# edge_graph(ax_26, gen26)#mc26)
edge_graph(ax_26, mc26)

ax_32 = Axis(fig[1, 2]; title="ρ=32", titlesize=36)
# edge_graph(ax_32, gen32)#mc32)
edge_graph(ax_32, mc32)


##
# 3D plot
# ax_new = LScene(fig[1, 2]; show_axis=false) # Axis3(fig[1,3])# 
# colors = [cgrad(colormap)[mc[i]] for i in eachindex(mc)]
# inds = 1:1:10000
# scatter!(ax_new, x[1, inds], x[2, inds], x[3, inds], color=colors[inds], markersize=20.0, markerspacing=0.1, markerstrokewidth=0.0)
# rotate_cam!(ax_new.scene, (0.0, -10.5, 0.0))

# save("figs/three-state-network-comparison-newest.png", fig)
display(fig)
