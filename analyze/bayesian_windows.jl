# generate plot of generator entries (mean holding times, actually) for sliding windows and KL div plot of the same

# here I need all of the static generators and all of the changing rho ones

# create dictionary of static environment markov chains
dt = 0.0 #inefficient?? all of this
static_mc_dict = Dict{Int64,Vector}()
for (i, ρₛ) in enumerate(range(26, 32)) 
    @info "opening data"
    hfile = h5open("data/lorenz"*string(ρₛ)*".hdf5")
    x = read(hfile["x"])
    dt = read(hfile["dt"]) #does this even work
    ρ = read(hfile["rho"])
    close(hfile)
    @info "applying embedding"
    markov_indices = zeros(Int, size(x,2))
    for i in ProgressBar(eachindex(markov_indices))
        markov_indices[i] = embedding(x[:, i])
    end
    static_mc_dict[ρ] = markov_indices
end

# create dict of sliding window chains?


#