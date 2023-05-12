# generate plot of mean holding times for sliding windows and KL div plot of the same
using Polynomials
using MarkovChainHammer.BayesianMatrix
using MarkovChainHammer.TransitionMatrix: steady_state, generator
using Distributions
using HDF5
using GLMakie
include("./analyze_util.jl")
# here I need all of the static sims and all of the changing rho ones

# create dictionary of static environment markov chains
dt = 0.0
static_mc_dict = Dict{Int64,Vector}()
for (i, ρₛ) in enumerate(range(26, 32)) 
    mc, dt = read_markov_chain("data/lorenz"*string(ρₛ)*".hdf5")
    static_mc_dict[ρₛ] = mc
end

# create list of sliding window chains
sliding_mc_dict = Dict()
for i in range(5,7)
    sliding_mcs = []
    mc, dt = read_markov_chain("./data/lorenz-changing-10e" * string(i) * ".hdf5")
    deltarho = Int64(floor(length(mc)/6))
    for j in range(0,4)
        push!(sliding_mcs, mc[deltarho*j+1:deltarho*(j+2)])
    end
    sliding_mc_dict[i] = sliding_mcs
end

middle_values = [27,28,29,30,31]
sliding_bayesian_dict = Dict()
for i in range(5,7)
    sliding_bayesians = [BayesianGenerator(mc; dt=dt) for mc in sliding_mc_dict[i]]
    sliding_bayesian_dict[i] = sliding_bayesians
end
sliding_reference = [BayesianGenerator(static_mc_dict[ρ]; dt=dt) for ρ in 26:32] #bayesian generators from static 

mc_delta, dt = read_markov_chain("./data/lorenz-changing-10e" * string(6) * ".hdf5")
bayesian_delta = BayesianGenerator(mc_delta; dt=dt)

#############################################

include("../visualize/window_plots.jl") # load in plotting functions themselves

function plot_diagonals(static_ref; sliding_windows=nothing, middle_values=nothing, bayesian_delta=nothing)
    ref_list = [x for x in 26:32]
    fig = Figure(resolution=(1200,800))
    ns = [1, 5, 9, 2, 6, 10]
    for i in 1:2, j in 1:3
        n = popfirst!(ns)
        ax = Axis(fig[i,j], title="$n")
        plot_evolution(n,n, ref_list, static_ref; sliding_windows=sliding_windows, middle_values=middle_values, bayesian_delta=bayesian_delta)
    end
    # for i in 1:3, j in 1:4
    #     box = j+4*(i-1) 
    #     ax = Axis(fig[i,j], title="$box")
    #     plot_evolution(box,box, ref_list, static_ref; sliding_windows=sliding_windows, middle_values=middle_values, bayesian_delta=bayesian_delta)
    # end
    fig
    # save("figs/diagonal_evolution.png", fig)
end

plot_diagonals(sliding_reference; sliding_windows=sliding_bayesian_dict[6], middle_values=middle_values, bayesian_delta=bayesian_delta)


##
function kl_div(p, q; significance=nothing)
    # for significance: integer indicating factor of sigma difference
    # p and q are passed in as two bayesian matrices
    m1, s1 = realize(p)
    m2, s2 = realize(q)
    metric = log.(s1./s2) .+ (s1.^2 .+ (m1.-m2).^2)./(2 .* s2.^2) .- 0.5
    if significance === nothing
        return metric
    else
        thresh =  log.(s1./s2) .+ (s1.^2 .+ (significance .* s1).^2)./(2 .* s2.^2) .- 0.5
        return metric .>= thresh
    end
end

function get_metrics(arg_list; significance=nothing) 
    lists = [[] for _ in 1:length(arg_list)]
    signif = [[] for _ in 1:length(arg_list)]
    u = 0
    for bayesian_list in arg_list #arg list ends up being a list of three lists of five bayseian generators
        u += 1
        for dist in bayesian_list[2:end]
            metric = kl_div(bayesian_list[1], dist) #first the reference, then one under question. also I'm not 100% clear on how it handles this
            push!(lists[u], metric)
            if significance !== nothing
                sig = kl_div(bayesian_list[1], dist; significance=significance)
                push!(signif[u], sig)
            end
        end
    end
    return lists, signif
end

metrics, metrics_sig = get_metrics([sliding_bayesian_dict[i] for i in 5:7]; significance=2)

function plot_diag_kl(metrics; metrics_sig)
    fig = Figure(resolution=(1200,800))
    ns = [1, 5, 9, 2, 6, 10]
    for i in 1:2, j in 1:3
        n = popfirst!(ns)
        ax = Axis(fig[i,j], title="$n")
        plot_kl(n, n, metrics; metrics_sig=metrics_sig)
    end
    # fig = Figure(resolution=(3200,1200))
    # for i in 1:3, j in 1:4
    #     box = j+4*(i-1) 
    #     ax = Axis(fig[i,j])
    #     plot_kl(box, box, metrics; metrics_sig=metrics_sig)
    #     axislegend(ax, position=:lt)
    # end
    # save("diagonal_kl.png", fig)
    fig
end
            
plot_diag_kl(metrics; metrics_sig)

## off-diagonal entries

function plot_off_diagonal(static_ref, metrics; metrics_sig=nothing, sliding_windows=nothing, middle_values=nothing, bayesian_delta=nothing)
    fig = Figure(resolution=(800,800))
    colors = ["red", "orange", "green", "blue", "violet"]
    ref_list = [x for x in 26:32]
    box_list = [(5,9), (8,12)]
    # plot evolution
    for n in eachindex(box_list)
        i, j = box_list[n]
        ax = Axis(fig[1,n], title="$i -> $j")
        plot_evolution(i, j, ref_list, static_ref; sliding_windows=sliding_windows, middle_values=middle_values, bayesian_delta=bayesian_delta)
    end
    # plot kl div
    for n in eachindex(box_list)
        i, j = box_list[n]
        ax = Axis(fig[2,n], title="$i -> $j")
        plot_kl(i,j,metrics;metrics_sig)
    end
    # save("figs/off_diagonal_evolution.png", fig)
    fig
end

plot_off_diagonal(sliding_reference, metrics;metrics_sig=metrics_sig, sliding_windows=sliding_bayesian_dict[6], middle_values=middle_values, bayesian_delta=bayesian_delta)


# a little experiment

steady_state(generator(static_mc_dict[27])) # example of the reference values
steady_state(generator(sliding_mc_dict[5][1]))

function discrete_kl(p,q)
    #p and q each a series of numbers, implied as a function of index
    return sum([p[i]*log(p[i]/q[i]) for i in eachindex(p)])
end

steady_states = [steady_state(generator(static_mc_dict[i])) for i in middle_values]

for i in eachindex(sliding_mc_dict[6])
    q = steady_state(generator(sliding_mc_dict[6][i]))
    kls = [discrete_kl(p,q) for p in steady_states]
    println(kls)
end

begin
    fig = Figure(resolution=(800,800))
    ax = Axis(fig[1,1])
    for i in eachindex(sliding_mc_dict[6])
        q = steady_state(generator(sliding_mc_dict[6][i]))
        kls = [discrete_kl(p,q) for p in steady_states]
        lines!(middle_values, kls, label="$(middle_values[i])")
    end
    axislegend(ax, position=:lt)
    fig
end