using HDF5, LorenzExtremeVisualization
using MarkovChainHammer.TransitionMatrix: holding_times

include("generate_util.jl")
include("../analyze/analyze_util.jl")

#generate holding times from extra long static reference runs: 10^8 timesteps (this is now the same as the references used)

for ρₛ in [26, 32]
    @info "generating data for static Lorenz ρ=$ρₛ"
    x, dt = lorenz_data(;timesteps=10^8, Δt = 0.005, res=1, ϵ=0.0, ρ=t -> ρₛ + 0 * t / (10^7 * 0.005))
    @info "applying embedding"
    markov_indices = zeros(Int, size(x,2))
    for i in ProgressBar(eachindex(markov_indices))
        markov_indices[i] = aggregated_embedding(x[:, i])
    end
    ht = holding_times(markov_indices, 3; dt=dt)
    ht = Vector{Vector{Float64}}(ht)
    @info "saving data"
    hfile = h5open(pwd() * "/data/lorenz" * string(ρₛ) * "_holding-times.hdf5", "w")
    hfile["dt"] = dt
    hfile["rho"] = ρₛ
    hfile["ht1"] = ht[1]
    hfile["ht2"] = ht[2]
    hfile["ht3"] = ht[3]
    close(hfile)
    @info "done saving holding times for static Lorenz ρ=$ρₛ"
end

#now from the the ensemble: 10^6 timesteps
begin
    @info "opening data"
    hfile = h5open(pwd() * "/data/lorenz-changing-10e6.hdf5")
    dt = read(hfile["dt"])
    mc = read(hfile["mc"])
    close(hfile)

    mc_ag = aggregate(mc)
    big_ht = holding_times(mc_ag, 3; dt=dt)

    hfile = h5open(pwd() * "/data/ensemble-members.hdf5")
    initials = read(hfile["initials"])
    ensemble_size = read(hfile["ensemble_size"])
    close(hfile)

    for i in 1:ensemble_size
        init = initials[:,i]
        @info "generating data for ensemble member $i"
        x, dt = lorenz_data(;timesteps=10^6, Δt = 0.005, res=1, ϵ=0.0, ρ=t -> 26 + 6 * t / (10^6 * 0.005), init=init)
        @info "applying embedding"
        markov_indices = zeros(Int, size(x,2))
        for j in ProgressBar(eachindex(markov_indices))
            markov_indices[j] = aggregated_embedding(x[:, j])
        end
        big_ht .= vcat.(big_ht, holding_times(markov_indices, 3; dt=dt)) #append holding times to big list 
    end

    big_ht = Vector{Vector{Float64}}(big_ht)

    @info "saving data"
    hfile = h5open(pwd() * "/data/ensemble-holding-times.hdf5", "w")
    hfile["ht1"] = big_ht[1]
    hfile["ht2"] = big_ht[2]
    hfile["ht3"] = big_ht[3]
    hfile["dt"] = dt
    close(hfile)
end


