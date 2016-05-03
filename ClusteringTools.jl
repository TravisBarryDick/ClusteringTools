module ClusteringTools

include("metrics.jl")
export FiniteMetric, PrecomputedMetric, DataMetric, precompute_metric, dist

include("clustering_optim.jl")
include("sparse_lp_solution.jl")
export SparseLPSolution, set_x!, get_x, set_y!, get_y, get_centers, get_cluster
export objective_value

include("round_ys.jl")
include("round_xs.jl")

include("convenience.jl")
export exact_clustering, rounded_clustering

include("distance_graphs.jl")
export DistanceGraph, add_edge!, node_group, total_size, make_metric

end
