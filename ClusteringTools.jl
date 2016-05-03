module ClusteringTools

include("utils.jl")

include("metrics.jl")
export FiniteMetric, PrecomputedMetric, DataMetric, precompute_metric, dist

include("cluster_sol.jl")
export ClusterSol, cluster_keys, num_clusters, num_points
export get_cluster, get_assignments, get_center, get_point, objective_value

include("kp_means.jl")
export kpmeans

include("clustering_optim.jl")
include("sparse_lp_solution.jl")
export SparseLPSolution, set_x!, get_x, set_y!, get_y


include("round_ys.jl")
include("round_xs.jl")

include("convenience.jl")
export exact_clustering, rounded_clustering

include("distance_graphs.jl")
export DistanceGraph, add_edge!, node_group, total_size, make_metric

end
