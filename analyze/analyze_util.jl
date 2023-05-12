function read_markov_chain(filename; include_x=false)
    @info "opening data"
    hfile = h5open(filename)
    dt = read(hfile["dt"])
    mc = read(hfile["mc"])
    if include_x
        x = read(hfile["x"])
    end
    close(hfile)
    if include_x
        return mc, dt, x
    end
    return mc, dt
end

function aggregate(markov_chain)
    new_mc = []
    for state in markov_chain
        if state <= 4
            push!(new_mc, 1)
        elseif state <= 8
            push!(new_mc, 2)
        else
            push!(new_mc, 3)
        end
    end
    return new_mc
end

function realize(Q)
    # Q is a bayesian matrix
    #return one matrix of means and one of stds; the diagonals are for holding times and off diagonals for exit probabilities
    number_of_states = size(Q)[1]
    rates = Q.posterior.rates
    exit_probabilities = Q.posterior.exit_probabilities
    real_means = zeros(number_of_states, number_of_states)
    real_stds = zeros(number_of_states, number_of_states)
    # handle diagonals first
    [real_means[i, i] = 1/mean(rates[i]) for i in eachindex(rates)]
    [real_stds[i, i] = real_means[i,i]^2*std(rates[i]) for i in eachindex(rates)] 
    # now off diagonals
    for i in 1:number_of_states
        alpha_j = exit_probabilities[i].alpha
        alpha_0 = sum(exit_probabilities[i].alpha)
        dist = NaN
        try
            dist = Distributions.Beta.(alpha_j, alpha_0 .- alpha_j)
        catch DomainError
            alpha_0 += 0.01 #THIS IS A HACK - THINK MORE !
            dist = Distributions.Beta.(alpha_j, alpha_0 .- alpha_j)
        end
        real_means[[1:i-1..., i+1:number_of_states...], i] = mean.(dist)
        real_stds[[1:i-1..., i+1:number_of_states...], i] = std.(dist)
    end
    return real_means, real_stds
end

