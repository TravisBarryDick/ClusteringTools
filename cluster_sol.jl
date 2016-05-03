type ClusterSol
    X::Array{Float64, 2}                 # The clustered data
    centers::Dict{Int, Vector{Float64}}  # Map from cluster keys to the centers
    p2cs::Vector{Set{Int}}               # Map from points to the cluster keys
    c2ps::Dict{Int, Set{Int}}            # Map from cluster keys to the points
end

"Returns a vector of the cluster keys."
cluster_keys(c::ClusterSol) = collect(keys(c.centers))

"Counts the number of clusters in the ClusterSol."
num_clusters(c::ClusterSol) = length(c.c2ps)

"Counts the number of points in the ClusterSol."
num_points(c::ClusterSol) = length(c.p2cs)

"Returns the points in cluster `key`."
get_cluster(c::ClusterSol, key) = c.c2ps[key]

"Returns the clusters keys that point `i` is assigned to."
get_assignments(c::ClusterSol, i) = c.p2cs[i]

"Returns the cluster with the given `key`."
get_center(c::ClusterSol, key) = c.centers[key]

"Returns the point at index `i`"
get_point(c::ClusterSol, i) = c.X[:,i]

_sqdist(x, y) = dot(x,x) + dot(y,y) - 2*dot(x,y)
_dist(x,y) = sqrt(_sqdist(x,y))

function objective_value(c::ClusterSol, obj)
    total = 0.0
    for i in 1:num_points(c)
        for key in get_assignments(c, i)
            if obj == :KMedian
                total += _dist(get_center(c, key), get_point(c, i))
            elseif obj == :KMeans
                total += _sqdist(get_center(c, key), get_point(c, i))
            end
        end
    end
    return total
end
