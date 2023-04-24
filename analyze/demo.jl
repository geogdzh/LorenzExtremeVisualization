# make holding time plot, generator, and network graph for static ρ=28 scenario 
include("lorenz_embedding.jl")

# here I need one static run, and then two static and one delta
markov_chain, dt = get_markov_chain("./data/lorenz28.hdf5")


## generator shown in paper:
Q_28 = generator(markov_chain; dt=dt)

## holding time histogram - cut from paper
# ht = holding_times(markov_chain)
# begin
#     fig = Figure(resolution=(1600, 1200))
#     for i in 1:3, j in 1:4
#         box = j+4*(i-1) 
#         ax = Axis(fig[i,j], title="$box") 
#         hist!(ax, ht[box], normalization=:probability,nbins=10)
#     end
#     fig
# end

## network graph

