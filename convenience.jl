function exact_clustering(data::Array{Float64,2}, k, p, ℓ, L, obj; kwargs...)
    metric = DataMetric(data)
    s = exact_clustering(metric, k, p, ℓ, L, obj; kwargs...)
    return ClusterSol(s, data)
end

function exact_clustering(metric::FiniteMetric, k, p, ℓ, L, obj; kwargs...)
    x, y = solve_clustering_optimization(metric, k, p, ℓ, L;
                                         obj=obj, solve_ip=true,
                                         kwargs...)
    return SparseLPSolution(x, y)
end

function rounded_clustering(data::Array{Float64,2}, k, p, ℓ, L, obj; kwargs...)
    metric = DataMetric(data)
    s = rounded_clustering(metric, k, p, ℓ, L, obj; kwargs...)
    return ClusterSol(s, data)
end

function rounded_clustering(metric::FiniteMetric, k, p, ℓ, L, obj; kwargs...)
    x, y = solve_clustering_optimization(metric, k, p, ℓ, L;
                                         obj=obj, solve_ip=false,
                                         kwargs...)
    s = SparseLPSolution(x, y)
    round_ys!(metric, s, obj)
    round_xs!(metric, s, p, ℓ, L)
    return s
end
