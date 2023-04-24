# generate plot of mean holding times for sliding windows and KL div plot of the same
using Polynomials
using Distributions
include("./lorenz_embedding.jl")
# here I need all of the static sims and all of the changing rho ones

# create dictionary of static environment markov chains
dt = 0.0 #inefficient?? all of this
static_mc_dict = Dict{Int64,Vector}()
for (i, ρₛ) in enumerate(range(26, 32)) 
    mc, dt = get_markov_chain("data/lorenz"*string(ρₛ)*".hdf5")
    static_mc_dict[ρₛ] = mc
end

# create list of sliding window chains
sliding_mc_dict = Dict()
for i in range(5,7)# range(5,7) # I don't think the 4 one is enough data at all
    sliding_mcs = []
    mc, dt = get_markov_chain("./data/lorenz-changing-10e" * string(i) * ".hdf5")
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

mc_delta, dt = get_markov_chain("./data/lorenz-changing-10e" * string(6) * ".hdf5")
bayesian_delta = BayesianGenerator(mc_delta; dt=dt)

##
function plot_evolution(n, k, ref_list, static_ref; sliding_windows=nothing, middle_values=nothing, bayesian_delta=nothing)
    #plot the static reference values of gen entries
    if n==k
        entry_list = [-1/(mean(static_ref[i])[n,k]) for i in eachindex(static_ref)] 
    else
        entry_list = Float64[]
        for z in eachindex(static_ref)
            gen = static_ref[z]
            alpha_j = gen.posterior.exit_probabilities[n].alpha[k-1] #j-1 only works if looking at transitions into later states!
            alpha_0 = sum(gen.posterior.exit_probabilities[n].alpha)
            dist = Distributions.Beta(alpha_j, alpha_0 - alpha_j)
            push!(entry_list,mean(dist))
        end 
    end
    scatter!(ref_list, entry_list)
    
    if sliding_windows !== nothing
        if n==k
            sliding_means = [-1/(mean(sliding_windows[i])[n,k]) for i in eachindex(sliding_windows)]
            sliding_stds = [(sliding_means[i]^2)*std(sliding_windows[i])[n,k] for i in eachindex(sliding_windows)] #according to error propagation
        else
            sliding_means = Float64[]
            sliding_stds = Float64[]
            for z in eachindex(sliding_windows)
                gen = sliding_windows[z]
                alpha_j = gen.posterior.exit_probabilities[n].alpha[k-1] #j-1 only works if looking at transitions into later states!
                alpha_0 = sum(gen.posterior.exit_probabilities[n].alpha)
                dist = Distributions.Beta(alpha_j, alpha_0 - alpha_j)
                push!(sliding_means,mean(dist))
                push!(sliding_stds, std(dist))
            end 
        end
        errorbars!(middle_values, sliding_means, sliding_stds, color="red")
        scatter!(middle_values, sliding_means, color="red")

        f1 = [Polynomials.fit(middle_values[1:2:end], sliding_means[1:2:end], 1)[i] for i in 0:1]
        f2 = [Polynomials.fit(middle_values[1:2:end], sliding_means[1:2:end], 2)[i] for i in 0:2]
        xs = ref_list[1]:0.1:ref_list[end]
        lines!(xs, [f1[1]+f1[2]*x for x in xs])
        lines!(xs, [f2[1]+f2[2]*x+f2[3]*x^2 for x in xs])
    end

    if bayesian_delta !== nothing
        if n==k
            deltamean = -1/(mean(bayesian_delta)[n,k]) 
            deltastd = (deltamean^2)*std(bayesian_delta)[n,k]
        else
            alpha_j = bayesian_delta.posterior.exit_probabilities[n].alpha[k-1] #j-1 only works if looking at transitions into later states!
            alpha_0 = sum(bayesian_delta.posterior.exit_probabilities[n].alpha)
            dist = Distributions.Beta(alpha_j, alpha_0 - alpha_j)
            deltamean = mean(dist)
            deltastd = std(dist)
        end
        scatter!(29, deltamean, color="blue")
        errorbars!([29], [deltamean], [deltastd], color="blue")
    end
end

function plot_diagonals(static_ref; sliding_windows=nothing, middle_values=nothing, bayesian_delta=nothing)
    ref_list = [x for x in 26:32]
    fig = Figure(resolution=(1600,1200))

    for i in 1:3, j in 1:4
        box = j+4*(i-1) 
        ax = Axis(fig[i,j], title="$box")
        plot_evolution(box,box, ref_list, static_ref; sliding_windows=sliding_windows, middle_values=middle_values, bayesian_delta=bayesian_delta)
    end
    fig
end

plot_diagonals(sliding_reference; sliding_windows=sliding_bayesian_dict[6], middle_values=middle_values, bayesian_delta=bayesian_delta)

##

function kl_div(p, q; significance=nothing)
    # for significance: integer indicating factor of sigma difference
    s1 = std(p)
    s2 = std(q)
    m1 = mean(p)
    m2 = mean(q)
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
    for bayesian_list in arg_list
        u += 1
        for dist in bayesian_list[2:end]
            metric = kl_div(bayesian_list[1], dist) #first the reference, then one under question
            push!(lists[u], metric)
            if significance !== nothing
                sig = kl_div(bayesian_list[1], dist; significance=significance)
                push!(signif[u], sig)
            end
        end
    end
    return lists, signif
end

# metrics, metrics_sig = get_metrics((sliding_bayesians, sl_bayesian_e5, sl_bayesian_e4); significance=2)
metrics, metrics_sig = get_metrics([sliding_bayesian_dict[i] for i in 5:7]; significance=2)

function plot_kl(n, k, metrics; metrics_sig)
    colors = ["red", "orange", "green", "blue", "violet"]
    labels = ["10e5", "10e6", "10e7"]
    xs = [x for x in 28:31]
    for list_no in eachindex(metrics) #iterate 1e6 through 1e4
        #all same color + label
        tmp = Float64[]
        marker_shapes = []
        for m in 1:4
            push!(tmp, metrics[list_no][m][n,k])
            if metrics_sig[list_no][m][n,k] 
                push!(marker_shapes, :circle) #o
            else
                push!(marker_shapes, :diamond) #s
            end
        end
        # marker_shapes = [tmp_sig[i] ? 'o' : 's' for i in 1:length(tmp_sig)]
        scatter!(xs, tmp, color=(colors[list_no], 0.5), marker=marker_shapes, markersize=20, strokewidth=2, strokecolor=:black)
        lines!(xs, tmp, color=colors[list_no], label=labels[list_no])
    end
end

function plot_diag_kl(metrics; metrics_sig)
    fig = Figure(resolution=(3200,1200))
    for i in 1:3, j in 1:4
        box = j+4*(i-1) 
        ax = Axis(fig[i,j])
        plot_kl(box, box, metrics; metrics_sig=metrics_sig)
        axislegend(ax, position=:lt)
    end
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
    fig
end

plot_off_diagonal(sliding_reference, metrics;metrics_sig=metrics_sig, sliding_windows=sliding_bayesian_dict[6], middle_values=middle_values, bayesian_delta=bayesian_delta)
