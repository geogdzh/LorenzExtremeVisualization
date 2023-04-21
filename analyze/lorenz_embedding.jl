using HDF5, GLMakie, LinearAlgebra, Statistics, Random, ProgressBars
using MarkovChainHammer.BayesianMatrix
using MarkovChainHammer.TransitionMatrix: steady_state, holding_times
using MarkovChainHammer.TransitionMatrix: perron_frobenius, generator
Random.seed!(12345)
##
# @info "opening data"
# hfile = h5open("data/lorenz28.hdf5")
# x = read(hfile["x"])
# dt = read(hfile["dt"])
# close(hfile)
##
function embedding(current_state)
    state = [0,0,0]
    state[1] = current_state[1] > 0 ? 1 : 0
    state[2] = current_state[2] > 0 ? 1 : 0
    state[3] = current_state[3] > 37.9 ? 2 : (current_state[3] > 22.5 ? 1 : 0)
    map = Dict([0,0,0] => 1, [1,0,0] => 2, [0,1,0] => 3, [1,1,0] => 4,
                [0,0,1] => 5, [1,0,1] => 6, [0,1,1] => 7, [1,1,1] => 8,
                [0,0,2] => 9, [1,0,2] => 10, [0,1,2] => 11, [1,1,2] => 12)
    return map[state]
end
##
# @info "applying embedding"
# markov_indices = zeros(Int, size(x,2))
# for i in ProgressBar(eachindex(markov_indices))
#     markov_indices[i] = embedding(x[:, i])
# end
##
# @info "generating Bayesian matrix"
# Q = mean(BayesianGenerator(markov_indices; dt = dt))

# gen26 = generator(markov_indices; dt=dt)
##
function get_markov_chain(filename)
    @info "opening data"
    hfile = h5open(filename)
    x = read(hfile["x"])
    dt = read(hfile["dt"])
    close(hfile)
    @info "applying embedding"
    markov_indices = zeros(Int, size(x,2))
    for i in ProgressBar(eachindex(markov_indices))
        markov_indices[i] = embedding(x[:, i])
    end
    return markov_indices, dt
end
