include("./generate_util.jl")

#generate static lorenz
for (i, ρₛ) in enumerate(range(26, 32)) 
    @info "generating data for static Lorenz ρ=$ρₛ"
    x, dt = lorenz_data(;timesteps=10^8, Δt = 0.005, res=1, ϵ=0.0, ρ=t -> ρₛ + 0 * t / (10^8 * 0.005)) #rerun!!
    @info "applying embedding"
    markov_indices = zeros(Int, size(x,2))
    for i in ProgressBar(eachindex(markov_indices))
        markov_indices[i] = embedding(x[:, i])
    end
    @info "saving data"
    hfile = h5open(pwd() * "/data/lorenz" * string(ρₛ) * ".hdf5", "w")
    hfile["x"] = x
    hfile["dt"] = dt
    hfile["rho"] = ρₛ
    hfile["mc"] = markov_indices
    close(hfile)
    @info "done saving data for static Lorenz ρ=$ρₛ"
end

# generate lorenz with time-varying rho
for i in range(5,7)
    @info "generating data for changing Lorenz at timesteps = 10^$i"
    x, dt = lorenz_data(;timesteps=10^i, Δt = 0.005, res=1, ϵ=0.0, ρ=t -> 26 + 6 * t / (10^i * 0.005)) 
    @info "applying embedding"
    markov_indices = zeros(Int, size(x,2))
    for i in ProgressBar(eachindex(markov_indices))
        markov_indices[i] = embedding(x[:, i])
    end
    @info "saving data"
    hfile = h5open(pwd() * "/data/lorenz-changing-10e" * string(i) * ".hdf5", "w")
    hfile["x"] = x
    hfile["dt"] = dt
    hfile["mc"] = markov_indices
    close(hfile)
    @info "done saving data for changing Lorenz at timesteps = 10^$i"
end

