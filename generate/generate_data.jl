using LorenzExtremeVisualization, HDF5, Random, ProgressBars
using MarkovChainHammer.TransitionMatrix: generator

# random seed for reproducibility
Random.seed!(12345)

# create data directory if it's not there
isdir(pwd() * "/data") ? nothing : mkdir(pwd() * "/data")

# generate Lorenz data
if isfile(pwd() * "/data/lorenz26.hdf5") #unideal because just checking for one
    @info "lorenz data already exists. skipping data generation"
else
    include("lorenz.jl")
end