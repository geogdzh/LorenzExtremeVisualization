using HDF5, GLMakie, LinearAlgebra, Statistics, Random
using MarkovChainHammer.BayesianMatrix
using MarkovChainHammer.TransitionMatrix: steady_state
using MarkovChainHammer.TransitionMatrix: perron_frobenius
using SparseArrays, Graphs, GraphMakie, Printf
##
include(pwd() * "/analyze/analyze_util.jl")

mc, dt, x = read_markov_chain("data/lorenz26.hdf5"; include_x=true)
mc26 = aggregate(mc)
mc, dt, x = read_markov_chain("data/lorenz32.hdf5"; include_x=true)
mc32 = aggregate(mc)

##

function edge_graph(ax, mc)
    Q = mean(BayesianGenerator(mc; dt = dt))
    # Generator
    g_Q = DiGraph(Q')
    Q_prim = zeros(size(Q))
    for i in 1:size(Q)[1]
        Q_prim[:, i] .= -Q[:, i] / Q[i, i]
        Q_prim[i, i] = -1 / Q[i, i]
    end

    elabels = string.([round(Q_prim[i]; digits=2) for i in 1:ne(g_Q)])
    # [@sprintf("%.0e", node_labels[i]) for i in 1:nv(G)]
    transparancy = [Q_prim[i] for i in 1:ne(g_Q)]
    elabels_color = [(:black, transparancy[i] > eps(100.0)) for i in 1:ne(g_Q)]
    # edge_color_Q = [(cgrad(colormap)[(i - 1)÷12 + 1], transparancy[i]) for i in 1:ne(g_Q)]
    edge_color_Q = [(colors[(i - 1)÷3 + 1], transparancy[i]) for i in 1:ne(g_Q)]
    # node_color = [(cgrad(colormap)[i]) for i in 1:nv(g_Q)]
    node_color = [(colors[i]) for i in 1:nv(g_Q)]
    edge_attr = (; linestyle=[:dot, :dash, :dash, :dash, :dot, :dash, :dash, :dash, :dot])
    elabels_fontsize = 18
    nlabels_fontsize = 18
    node_size = 60.0

    edge_width_Q = [5.0 for i in 1:ne(g_Q)]
    arrow_size_Q = [20.0 for i in 1:ne(g_Q)]
    node_labels_Q = repr.(1:nv(g_Q))

    #                   ; edge_color=edge_color_Q, edge_width=edge_width_Q
    kwargs_edges = (; elabels=elabels, elabels_color=elabels_color, elabels_fontsize=elabels_fontsize, edge_color=edge_color_Q, edge_width=edge_width_Q)
    kwargs_nodes = (; node_color=node_color, node_size=node_size, nlabels=node_labels_Q, nlabels_fontsize=nlabels_fontsize)
    kwargs_arrows = (; arrow_size=arrow_size_Q)
        
    p_Q = graphplot!(ax, g_Q; kwargs_edges..., kwargs_nodes..., kwargs_arrows...)

    hidedecorations!(ax)
end


fig = Figure(resolution=(2000, 1000))
colors = ["blue", "violet", "red"]

ax_26 = Axis(fig[1, 1]; title="ρ=26", titlesize=30)
edge_graph(ax_26, mc26)

ax_32 = Axis(fig[1, 2]; title="ρ=32", titlesize=30)
edge_graph(ax_32, mc32)


##
# 3D plot
# ax_new = LScene(fig[1, 2]; show_axis=false) # Axis3(fig[1,3])# 
# colors = [cgrad(colormap)[mc[i]] for i in eachindex(mc)]
# inds = 1:1:10000
# scatter!(ax_new, x[1, inds], x[2, inds], x[3, inds], color=colors[inds], markersize=20.0, markerspacing=0.1, markerstrokewidth=0.0)
# rotate_cam!(ax_new.scene, (0.0, -10.5, 0.0))

save("figs/three-state-network-comparison.png", fig)
display(fig)
