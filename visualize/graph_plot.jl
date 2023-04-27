using HDF5, GLMakie, LinearAlgebra, Statistics, Random
using MarkovChainHammer.BayesianMatrix
using MarkovChainHammer.TransitionMatrix: steady_state
using MarkovChainHammer.TransitionMatrix: perron_frobenius
using SparseArrays, Graphs, GraphMakie, Printf
##
include(pwd() * "/analyze/lorenz_embedding.jl")
ρₛ = 26
@info "opening data"
file = "lorenz" * string(ρₛ) * ".hdf5"
hfile = h5open("data/" * file )
x = read(hfile["x"])
dt = read(hfile["dt"])
X, dt = get_markov_chain("data/lorenz" * string(ρₛ) * ".hdf5")
close(hfile)
##
Q = mean(BayesianGenerator(X; dt = dt))
fig = Figure(resolution=(2000, 1000))
colormap = :glasbey_hv_n256
ax_Q = Axis(fig[1, 1]; title="Generator", titlesize=30)

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
edge_color_Q = [(:black, transparancy[i]) for i in 1:ne(g_Q)]
node_color = [(cgrad(colormap)[i]) for i in 1:nv(g_Q)]
edge_attr = (; linestyle=[:dot, :dash, :dash, :dash, :dot, :dash, :dash, :dash, :dot])
elabels_fontsize = 18
nlabels_fontsize = 18
node_size = 30.0

edge_width_Q = [5.0 for i in 1:ne(g_Q)]
arrow_size_Q = [20.0 for i in 1:ne(g_Q)]
node_labels_Q = repr.(1:nv(g_Q))

kwargs_edges = (; elabels=elabels, elabels_color=elabels_color, elabels_fontsize=elabels_fontsize, edge_color=edge_color_Q, edge_width=edge_width_Q)
kwargs_nodes = (; node_color=node_color, node_size=node_size, nlabels=node_labels_Q, nlabels_fontsize=nlabels_fontsize)
kwargs_arrows = (; arrow_size=arrow_size_Q)
     
p_Q = graphplot!(ax_Q, g_Q; kwargs_edges..., kwargs_nodes..., kwargs_arrows...)

hidedecorations!(ax_Q)

##
# 3D plot
ax_new = LScene(fig[1, 2]; show_axis=false) # Axis3(fig[1,3])# 
colors = [cgrad(colormap)[X[i]] for i in eachindex(X)]
inds = 1:1:10000
scatter!(ax_new, x[1, inds], x[2, inds], x[3, inds], color=colors[inds], markersize=20.0, markerspacing=0.1, markerstrokewidth=0.0)
rotate_cam!(ax_new.scene, (0.0, -3.5, 0.0))
display(fig)