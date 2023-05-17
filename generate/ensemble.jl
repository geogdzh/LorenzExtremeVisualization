using LorenzExtremeVisualization, HDF5, Random, ProgressBars
using MarkovChainHammer.TransitionMatrix: generator

#### generate ensemble initial conditions
ensemble_size = 100

initials = zeros(3, ensemble_size)

for i in 1:ensemble_size
    init = [14.0, 15.0, 27.0] .+ rand(Float64, (3))
    x, dt = lorenz_data(;timesteps=10^5, Δt=0.005, res=10, ϵ=0.0, ρ=t -> 26.0, init=init)
    initials[:,i] = x[:, end]
end

hfile = h5open(pwd() * "/data/ensemble-members.hdf5", "w")
hfile["initials"] = initials
hfile["ensemble_size"] = ensemble_size
close(hfile)