include("./lorenz_embedding.jl")

function lorenz(x, ρ, σ, β)
    x1 = x[1]
    x2 = x[2]
    x3 = x[3]
    return [σ * (x2 - x1), x1 * (ρ - x3) - x2, x1 * x2 - β * x3]
end

function lorenz_data(;timesteps=10^7, Δt=0.005, res=10, ϵ=0.0, ρ=28.0, init=[14.0, 15.0, 27.0]) # change default for rho
    rhs(x, t) = lorenz(x, ρ(t), 10.0, 8.0 / 3.0)
    x_f = zeros(3, timesteps)
    x_f[:, 1] = init#[14.0, 15.0, 27.0]
    evolve = RungeKutta4(3)
    for i in ProgressBar(2:timesteps)
        xOld = x_f[:, i-1]
        evolve(rhs, xOld, Δt)
        𝒩 = randn(3)
        @inbounds @. x_f[:, i] = evolve.xⁿ⁺¹ + ϵ * sqrt(Δt) * 𝒩
    end
    L2 = floor(Int, timesteps / res)
    Dt = Δt * res
    x = zeros(3, L2)
    for i in 1:L2
        @inbounds x[:, i] .= x_f[:, res*i]
    end

    return x, Dt
end

# lorenz_data(; timesteps=10^7, Δt=0.005, res=1, ϵ=0.0, ρ = t -> 28.0 + 6 * t / (timesteps * Δt)) = lorenz_data(timesteps, Δt, res, ϵ, ρ)

for (i, ρₛ) in enumerate(range(26, 32)) 
    @info "generating data for static Lorenz ρ=$ρₛ"
    x, dt = lorenz_data(;timesteps=10^7, Δt = 0.005, res=1, ϵ=0.0, ρ=t -> ρₛ + 0 * t / (10^7 * 0.005))
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


begin
    # @info "opening data"
    hfile = h5open(pwd() * "/data/lorenz-changing-10e6.hdf5")
    dt = read(hfile["dt"])
    mc = read(hfile["mc"])
    close(hfile)
    generators = generator(mc; dt=dt)

    Q_slices = [] #not the best way to do it
    deltarho = Int64(floor(length(mc)/6))
    for j in range(0,4)
        mc_slice = mc[deltarho*j+1:deltarho*(j+2)]
        gen = generator(mc_slice; dt=dt)
        push!(Q_slices, gen)
    end

    hfile = h5open(pwd() * "/data/ensemble-members.hdf5")
    intials = read(hfile["initials"])
    ensemble_size = read(hfile["ensemble_size"])
    close(hfile)

    for i in 1:ensemble_size
        init = initials[:,i]
        # init = [14.0, 15.0, 27.0] .+ (i * 0.25 * (-1)^i)
        @info "generating data for ensemble member $i"
        x, dt = lorenz_data(;timesteps=10^6, Δt = 0.005, res=1, ϵ=0.0, ρ=t -> 26 + 6 * t / (10^6 * 0.005), init=init)
        @info "applying embedding"
        markov_indices = zeros(Int, size(x,2))
        for j in ProgressBar(eachindex(markov_indices))
            markov_indices[j] = embedding(x[:, j])
        end
        gen = generator(markov_indices; dt=dt)
        generators = [generators gen]
        @info "splitting markov chain into windows"
        deltarho = Int64(floor(length(markov_indices)/6))
        for j in range(0,4)
            mc_slice = markov_indices[deltarho*j+1:deltarho*(j+2)]
            gen = generator(mc_slice; dt=dt)
            Q_slices[j+1] = [Q_slices[j+1] gen]
        end
    end
    @info "saving data"
    hfile = h5open(pwd() * "/data/lorenz-changing-ensemble.hdf5", "w")
    hfile["Q_full"] = generators
    for i in 1:5 #also feels unideal?
        hfile["Q$i"] = Q_slices[i]
    end
    hfile["dt"] = dt
    hfile["ensemble_size"] = ensemble_size
    close(hfile)
end