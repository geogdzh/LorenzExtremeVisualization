# generate all plots deriving from sliding-window timeseries grouping
# using Polynomials
using MarkovChainHammer.BayesianMatrix
using MarkovChainHammer.TransitionMatrix: steady_state, generator
using Distributions
using HDF5
using GLMakie
include("./analyze_util.jl")
include("../generate/generate_util.jl")

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

#create dict of sliding window generators
middle_values = [27,28,29,30,31]
sliding_bayesian_dict = Dict()
for i in range(5,7)
    sliding_bayesians = [BayesianGenerator(mc; dt=dt) for mc in sliding_mc_dict[i]]
    sliding_bayesian_dict[i] = sliding_bayesians
end
sliding_reference = [BayesianGenerator(static_mc_dict[ρ]; dt=dt) for ρ in 26:32] #bayesian generators from static 

# singular overall changing average
mc_delta, dt = read_markov_chain("./data/lorenz-changing-10e" * string(6) * ".hdf5")
bayesian_delta = BayesianGenerator(mc_delta; dt=dt)

# load in ensemble
ensemble_generators = [] 
for i in 1:5
    gen = load_Q_bayes("./data/ensemble-generators/Q$i.hdf5")
    push!(ensemble_generators, gen)
end
push!(ensemble_generators, load_Q_bayes("./data/ensemble-generators/Q_full.hdf5"))

hfile = h5open("./data/ensemble-generators/Q_full.hdf5")
ensemble_size = read(hfile["ensemble_size"])
close(hfile)


#############################################

include("../visualize/window_plots.jl") # load in plotting functions themselves

function plot_diagonals(static_ref; sliding_windows=nothing, middle_values=nothing, bayesian_delta=nothing, ensemble_generators=nothing)
    ref_list = [x for x in 26:32]
    fig = Figure(resolution=(1200,800))
    ns = [1, 5, 9, 2, 6, 10]
    for i in 1:2, j in 1:3
        n = popfirst!(ns)
        ax = Axis(fig[i,j], title="State $n", titlesize=22, ylabel=(j==1 ? "Mean holding time" : ""), ylabelsize=20, xlabel=(i==2 ? "ρ" : ""), xlabelsize=20)#, xticklabelsize=18, yticklabelsize=18)
        plot_evolution(n,n, ref_list, static_ref; sliding_windows=sliding_windows, middle_values=middle_values, bayesian_delta=bayesian_delta, ensemble_generators=ensemble_generators)
        if i == 1
            hidexdecorations!(ax, grid=false)
        end
        if n == 1
            axislegend(ax)
        end

    end
    ### select below instead for all 12 states
    # for i in 1:3, j in 1:4
    #     box = j+4*(i-1) 
    #     ax = Axis(fig[i,j], title="$box")
    #     plot_evolution(box,box, ref_list, static_ref; sliding_windows=sliding_windows, middle_values=middle_values, bayesian_delta=bayesian_delta)
    # end
    fig
    save("figs/diagonal_evolution_ensemble.png", fig)
end

plot_diagonals(sliding_reference; sliding_windows=sliding_bayesian_dict[6], middle_values=middle_values, ensemble_generators=ensemble_generators) #bayesian_delta=bayesian_delta

# gives errors of 0.029391288674768482 for linear and 0.027833771539007426 for quadratic

##
function kl_div(p, q; significance=nothing)
    # significance: integer indicating factor of sigma difference
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
            metric = kl_div(bayesian_list[1], dist)
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
    xticks = ["(27,29)", "(28,30)", "(29,31)", "(30,32)"]
    for i in 1:2, j in 1:3
        n = popfirst!(ns)
        ax = Axis(fig[i,j], title="State $n", ylabel=(j==1 ? "KL-div" : ""), ylabelsize=20, xlabel=(i==2 ? "ρ" : ""), xlabelsize=20, xticks=(28:31, xticks))
        plot_kl(n, n, metrics; metrics_sig=metrics_sig)
        if i == 1
            hidexdecorations!(ax, grid=false)
        end
        if n == 1
            axislegend(ax, position=:lt)
        end
    end
    # fig = Figure(resolution=(3200,1200))
    # for i in 1:3, j in 1:4
    #     box = j+4*(i-1) 
    #     ax = Axis(fig[i,j])
    #     plot_kl(box, box, metrics; metrics_sig=metrics_sig)
    #     axislegend(ax, position=:lt)
    # end
    save("figs/diagonal_kl.png", fig)
    fig
end
            
plot_diag_kl(metrics; metrics_sig)

## off-diagonal entries

function plot_off_diagonal(static_ref, metrics; metrics_sig=nothing, sliding_windows=nothing, middle_values=nothing, bayesian_delta=nothing, ensemble_generators=nothing)
    fig = Figure(resolution=(800,800))
    colors = ["red", "orange", "green", "blue", "violet"]
    ref_list = [x for x in 26:32]
    box_list = [(5,9), (8,12)]
    xticks = ["(27,29)", "(28,30)", "(29,31)", "(30,32)"]
    # plot evolution
    for n in eachindex(box_list)
        i, j = box_list[n]
        ax = Axis(fig[1,n], title="Transition $i -> $j",  ylabel=(n==1 ? "Exit probability" : ""), ylabelsize=20, titlesize=20)
        plot_evolution(i, j, ref_list, static_ref; sliding_windows=sliding_windows, middle_values=middle_values, bayesian_delta=bayesian_delta, ensemble_generators=ensemble_generators)
        if n==2
            axislegend(ax, position=:rb)
        end
    end
    # plot kl div
    for n in eachindex(box_list)
        i, j = box_list[n]
        ax = Axis(fig[2,n], title="Transition $i -> $j", ylabel=(n==1 ? "KL-div" : ""), ylabelsize=20, xlabel="ρ", xlabelsize=20, titlesize=20, xticks=(28:31, xticks))
        plot_kl(i,j,metrics;metrics_sig)
        if n==1
            axislegend(ax, position=:lt)
        end
    end
    save("figs/newest/off_diagonal_evolution_ensemble-100.png", fig)
    fig
end

plot_off_diagonal(sliding_reference, metrics;metrics_sig=metrics_sig, sliding_windows=sliding_bayesian_dict[6], middle_values=middle_values,  ensemble_generators=ensemble_generators) #bayesian_delta=bayesian_delta,

# gives errors of avg. 0.01807270209155247 for linear and 0.016877226649045382 for quadratic

####### Steady state
function discrete_kl(p,q)
    #p and q each a series of numbers, implied as a function of index
    return sum([p[i]*log(p[i]/q[i]) for i in eachindex(p)])
end


full_values = [26, middle_values...,32]
steady_states = [steady_state(generator(static_mc_dict[i])) for i in full_values]

colors = [:red, :orange, :lightblue, :blue, :violet]
begin
    fig = Figure(resolution=(800,600))
    ax = Axis(fig[1,1], xlabel="ρ (reference)", ylabel="KL-div")#, ylabelsize=26, xlabelsize=26, xticklabelsize=20, yticklabelsize=20)
    for i in eachindex(sliding_mc_dict[6]) #it's just the five windows
        q = steady_state(sliding_bayesian_dict[6][i]) # this is for the original single-instance one
        kls = [discrete_kl(p,q) for p in steady_states]
        lines!(full_values, kls, label="$(middle_values[i]-1) - $(middle_values[i]+1)", color=colors[i])
        scatter!(full_values, kls, color=colors[i])

        for j in 2:ensemble_size # essentially need n generators for whatever the size of the ensemble is
            q = steady_state(rand(ensemble_generators[i]))
            kls = [discrete_kl(p,q) for p in steady_states]
            lines!(full_values, kls, color=(colors[i],0.1), linestyle=:dot)
        end
    end
    axislegend("Average over ρ=", position=:ct)#, labelsize=20, titlesize=26)
    save("figs/steady_state_kl_big.png", fig)
    fig
end