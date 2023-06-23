function get_means_stds(list_of_bayesians, n, k)
    if n==k
        means = [-1/(mean(list_of_bayesians[i])[n,k]) for i in eachindex(list_of_bayesians)]
        stds = [(means[i]^2)*std(list_of_bayesians[i])[n,k] for i in eachindex(list_of_bayesians)] #according to error propagation
    else
        means = Float64[]
        stds = Float64[]
        for z in eachindex(list_of_bayesians)
            gen = list_of_bayesians[z]
            alpha_j = gen.posterior.exit_probabilities[n].alpha[k-1] #j-1 only works if looking at transitions into later states!
            alpha_0 = sum(gen.posterior.exit_probabilities[n].alpha)
            dist = Distributions.Beta(alpha_j, alpha_0 - alpha_j)
            push!(means,mean(dist))
            push!(stds, std(dist)) # this is all good bc we do just want to look at the exit probabilities
        end 
    end
    return means, stds
end



function plot_evolution(n, k, ref_list, static_ref; sliding_windows=nothing, middle_values=nothing, bayesian_delta=nothing, ensemble_generators=nothing)
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
    scatter!(ref_list, entry_list, color=:black)
    
    if sliding_windows !== nothing
        sliding_means, sliding_stds = get_means_stds(sliding_windows, n, k)
        errorbars!(middle_values, sliding_means, sliding_stds.*2, color="red")
        shapes = [:star5, :circle, :star5, :circle, :star5]
        scatter!(middle_values, sliding_means, color="red", marker=shapes)
        # ** need to fix this to have a normalization if it's off diagonal !!!!
        f1 = [Polynomials.fit(middle_values[1:2:end], sliding_means[1:2:end], 1)[i] for i in 0:1]
        f2 = [Polynomials.fit(middle_values[1:2:end], sliding_means[1:2:end], 2)[i] for i in 0:2]
        xs = ref_list[1]:0.1:ref_list[end]
        lines!(xs, [f1[1]+f1[2]*x for x in xs])
        lines!(xs, [f2[1]+f2[2]*x+f2[3]*x^2 for x in xs])
        errors = []
        for pair in [(2,28),(4,30)]
            i, x = pair
            err1 = abs((f1[1]+f1[2]*x)-sliding_means[i])/sliding_means[i]
            err2 = abs((f2[1]+f2[2]*x+f2[3]*x^2)-sliding_means[i])/sliding_means[i]
            push!(errors, err1)
            push!(errors, err2)
        end
    end

    if bayesian_delta !== nothing #REFACTOR THIS with get_means_stds
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
        errorbars!([29], [deltamean], [deltastd*2], color="blue")
    end

    if !isnothing(ensemble_generators)
        ensemble_means, ensemble_stds = get_means_stds(ensemble_generators, n, k)
        ensemble_means = ensemble_means[1:5]
        ensemble_stds = ensemble_stds[1:5]
        band!(middle_values, ensemble_means .- ensemble_stds .* 2, ensemble_means .+ ensemble_stds .* 2, color=(:red, 0.3) )
    end

    @info "printing errors for $n,$k"
    println(errors)
    # return errors
end


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


