include("./generate_util.jl")

# generate ensemble
begin
    @info "opening data"
    hfile = h5open(pwd() * "/data/lorenz-changing-10e6.hdf5")
    dt = read(hfile["dt"])
    mc = read(hfile["mc"])
    close(hfile)

    current_generators = []

    deltarho = Int64(floor(length(mc)/6))
    for j in range(0,4)
        mc_slice = mc[deltarho*j+1:deltarho*(j+2)]
        gen = BayesianGenerator(mc_slice; dt=dt)
        push!(current_generators, gen)
    end
    push!(current_generators, BayesianGenerator(mc;dt=dt))

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
            markov_indices[j] = embedding(x[:, j])
        end
        
        @info "splitting markov chain into windows and updating generators"
        deltarho = Int64(floor(length(markov_indices)/6))
        for j in range(0,4)
            mc_slice = markov_indices[deltarho*j+1:deltarho*(j+2)]
            current_generators[j+1] = BayesianGenerator(mc_slice, current_generators[j+1].posterior; dt=dt)
        end
        current_generators[6] = BayesianGenerator(markov_indices, current_generators[6].posterior; dt=dt)
    end

    @info "saving data"
    for i in 1:5 
        save_Q_bayes(current_generators[i], "./data/ensemble-generators/Q$i.hdf5"; ensemble_size=ensemble_size, dt=dt) 
    end
    save_Q_bayes(current_generators[6], "./data/ensemble-generators/Q_full.hdf5"; ensemble_size=ensemble_size, dt=dt) 
end