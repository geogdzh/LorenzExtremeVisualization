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

# generate ensemble initial conditions
if isfile(pwd() * "/data/ensemble-members.hdf5") 
    @info "ensemble member list already exists. skipping initial condition generation"
else
    include("ensemble.jl")
end

# generate ensemble data
isdir(pwd() * "/data/ensemble-generators") ? nothing : mkdir(pwd() * "/data/ensemble-generators")

if isfile(pwd() * "/data/ensemble-generators/Q_full.hdf5") 
    @info "lorenz ensemble data already exists. skipping data generation"
else
    include("ensemble_lorenz.jl")
end

# generate holding times from extra-long simulations
if isfile(pwd() * "/data/lorenz26_holding-times.hdf5") 
    @info "extra holding times data already exists. skipping data generation"
else
    include("extra_holding_times.jl")
end