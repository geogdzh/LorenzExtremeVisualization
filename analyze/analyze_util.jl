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