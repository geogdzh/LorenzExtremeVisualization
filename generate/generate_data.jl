using LorenzExtremeVisualization, HDF5, Random, ProgressBars

# random seed for reproducibility
Random.seed!(12345)

# create data directory if it's not there
isdir(pwd() * "/data") ? nothing : mkdir(pwd() * "/data")

# potential well
if isfile(pwd() * "/data/lorenz.hdf5")
    @info "lorenz data already exists. skipping data generation"
else
    include("lorenz.jl")
end