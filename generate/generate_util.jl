using HDF5
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

#####
#ideally these should be recursive...

function save_Q_bayes(Q_bayes, filename; ensemble_size=nothing, dt=nothing)
    parameters = unpack(Q_bayes)
    hfile = h5open(filename, "w")
    outer_group = create_group(hfile, "parameters")
    for dist in keys(parameters)
        dist_group = create_group(outer_group, String(dist)) #iterating over dist = prior, posterior, predictive
        nested = getfield(parameters, dist)
        for field in keys(nested)
            value = getfield(nested, field)
            dist_group[String(field)] = value
        end
    end
    if !isnothing(ensemble_size)
        hfile["ensemble_size"] = ensemble_size
        hfile["dt"] = dt
    end
    close(hfile)
end


function load_Q_bayes(filename; has_ensemble_size=false, has_dt=false)
    hfile = h5open(filename)
    outer_names = [] #prior, posterior, predictive
    outer_values = []
    for dist in keys(hfile["parameters"]) #prior, posterior, predictive
        names = [] 
        values = []
        for key in keys(hfile["parameters"][dist])
            push!(names, Symbol(key))
            push!(values, read(hfile["parameters"][dist][key]))
        end
        innertup = (; zip(names, values)...)
        push!(outer_names, Symbol(dist))
        push!(outer_values, innertup)
    end
    parameters = (; zip(outer_names, outer_values)...)
    
    if has_ensemble_size #CHANGE so these are separate; could this be improved w dispatch?
        ensemble_size = read(hfile["ensemble_size"])
        dt = read(hfile["dt"])
        close(hfile)
        return BayesianGenerator(parameters), ensemble_size, dt
    end
    close(hfile)
    return BayesianGenerator(parameters)
end