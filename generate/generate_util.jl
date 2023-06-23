using HDF5

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