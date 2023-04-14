using HDF5, GLMakie, LinearAlgebra, Statistics, Random
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
lines(x)