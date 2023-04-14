using HDF5, GLMakie, LinearAlgebra, Statistics, Random, ProgressBars
using MarkovChainHammer.BayesianMatrix
using MarkovChainHammer.TransitionMatrix: steady_state
using MarkovChainHammer.TransitionMatrix: perron_frobenius
Random.seed!(12345)
##
@info "opening data"
hfile = h5open("data/lorenz.hdf5")
x = read(hfile["x"])
dt = read(hfile["dt"])
close(hfile)
##
# This is where the partitioning strategy should be implemented (won't need markov states for your method)
# This should just be a mapping for R³ -> {1, 2, ..., N} where N is the number of states
function embedding(current_state)
    markov_states = [[-sqrt(72), -sqrt(72), 27], [0.0, 0.0, 0.0], [sqrt(72), sqrt(72), 27]]
    argmin([norm(current_state - markov_state) for markov_state in markov_states])
end
##
@info "applying embedding"
markov_indices = zeros(Int, size(x,2))
for i in ProgressBar(eachindex(markov_indices))
    markov_indices[i] = embedding(x[:, i])
end
##
@info "generating Bayesian matrix"
Q = mean(BayesianGenerator(markov_indices; dt = dt))