function lorenz(x, ρ, σ, β)
    x1 = x[1]
    x2 = x[2]
    x3 = x[3]
    return [σ * (x2 - x1), x1 * (ρ - x3) - x2, x1 * x2 - β * x3]
end

function lorenz_data(;timesteps=10^7, Δt=0.005, res=10, ϵ=0.0, ρ=28.0)
    rhs(x, t) = lorenz(x, ρ(t), 10.0, 8.0 / 3.0)
    x_f = zeros(3, timesteps)
    x_f[:, 1] = [14.0, 15.0, 27.0]
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
    x, dt = lorenz_data(;timesteps=10^7, Δt = 0.005, res=10, ϵ=0.0, ρ=t -> ρₛ + 0 * t / (10^7 * 0.005))
    @info "saving data for static Lorenz ρ=$ρₛ"
    hfile = h5open(pwd() * "/data/lorenz" * string(ρₛ) * ".hdf5", "w")
    hfile["x"] = x
    hfile["dt"] = dt
    hfile["rho"] = ρₛ
    close(hfile)
    @info "done saving data for static Lorenz ρ=$ρₛ"
end

for i in range(4,7) #exact range tbd
    x, dt = lorenz_data(;timesteps=10^i, Δt = 0.005, res=10, ϵ=0.0, ρ=t -> 26 + 6 * t / (10^i * 0.005)) 
    @info "saving data for changing Lorenz at timesteps = 10^$i"
    hfile = h5open(pwd() * "/data/lorenz-changing-10e" * string(i) * ".hdf5", "w")
    hfile["x"] = x
    hfile["dt"] = dt
    hfile["rho"] = 26
    close(hfile)
    @info "done saving data for changing Lorenz at timesteps = 10^$i"
end