import Clustering

# An implementation of k-means++ with p replication

"Represents the state during the Lloyd iteration of k-means++"
type LloydState
    centers::Array{Float64, 2} # Matrix whose columns are centers
    assignments::Array{Int, 2} # Matrix whose columns are assignment vectors
    counts::Array{Int, 1}      # Counts of how many points belong to each center
end

"Set the center of each cluster to be the center of mass of its points."
function update_centers!(X, s::LloydState)
    d, n = size(X)
    p = size(s.assignments, 1)
    fill!(s.centers, 0.0)
    for i in 1:n
        for j in 1:p
            add_cols!(s.centers, s.assignments[j,i], X, i)
        end
    end
    for i in 1:size(s.centers,2)
        scale_col!(s.centers, i, 1.0 / s.counts[i])
    end
end

"Assign each point to its p closest centers."
function update_assignments!(X, s::LloydState)
    d, n = size(X)
    p = size(s.assignments, 1)
    k = size(s.centers, 2)
    fill!(s.counts, 0)
    sqdists = Array(Float64, k)
    any_changed = false
    for i in 1:n
        for j in 1:k
            sqdists[j] = col_sqdist(X, i, s.centers, j)
        end
        pth_smallest_distance = select(sqdists, p)
        cix_count = 1
        for j in 1:k
            if sqdists[j] <= pth_smallest_distance
                any_changed |= (s.assignments[cix_count, i] != j)
                s.assignments[cix_count, i] = j
                s.counts[j] += 1
                cix_count += 1
            end
            if cix_count > p
                break
            end
        end
    end
    any_changed
end

"Choose initial centers according to k-means++."
function kmpp_init(X, k)
    n = size(X,2)
    centers = Array(Int, k)
    centers[1] = rand(1:n)
    sqdists = Array(Float64, n)
    total_sqdist = 0.0
    for i in 1:n
        sqdists[i] = col_sqdist(X, i, X, centers[1])
        total_sqdist += sqdists[i]
    end
    for i in 2:k
        centers[i] = discrete_sample(sqdists, total_sqdist)
        for j in 1:n
            sqdist = col_sqdist(X, j, X, centers[i])
            if sqdist < sqdists[j]
                total_sqdist += sqdist - sqdists[j]
                sqdists[j] = sqdist
            end
        end
    end
    return centers
end

"""
Clusters the data using k-means++ initialization and a modified Lloyd's method
that maintains p assignments for each point.
"""
function kpmeans(X, k, p)
    d, n = size(X)
    # Get the starting centers from kmeans++ initialization
    cixs = kmpp_init(X, k)
    centers = Array(Float64, d, k)
    for i in 1:k
        copy_col!(centers, i, X, cixs[i])
    end
    # Construct the solution with the given centers and 0 assignments
    sol = LloydState(centers, zeros(Int,p,n), zeros(Int,k))
    # Do Lloyd's iteration until assignments converge
    while update_assignments!(X, sol)
        update_centers!(X, sol)
    end
    return ClusterSol(sol, X)
end

"Constructs a `ClusterSol` from a `LloydState`."
function ClusterSol(s::LloydState, X::Array{Float64, 2})
    (p,n) = size(s.assignments)
    k = size(s.centers, 2)
    centers = [i => s.centers[:,i] for i in 1:k]
    p2cs = [Set{Int}() for i in 1:n]
    c2ps = [i => Set{Int}() for i in 1:k]
    for i in 1:n
        for j in 1:p
            push!(p2cs[i], s.assignments[j,i])
            push!(c2ps[s.assignments[j,i]], i)
        end
    end
    return ClusterSol(X, centers, p2cs, c2ps)
end

"Constructs a `ClsuterSol` from a `Clsutering.KmeansResult`."
function ClusterSol(s::Clustering.KmeansResult, X::Array{Float64, 2})
    n = length(s.assignments)
    k = size(s.centers,2)
    centers = [i => s.centers[:,i] for i in 1:k]
    p2cs = [Set{Int}() for i in 1:n]
    c2ps = [i => Set{Int}() for i in 1:k]
    for i in 1:n
        push!(p2cs[i], s.assignments[i])
        push!(c2ps[s.assignments[i]], i)
    end
    return ClusterSol(X, centers, p2cs, c2ps)
end

"Sets a cluster's center to be the center of mass of the points assined to it."
function update_center!(c::ClusterSol, key)
    d = size(c.X,1)
    c.centers[key] = zeros(d)
    for p in c.c2ps[key]
        c.centers[key] += slice(c.X, :, p)
    end
    c.centers[key] /= length(c.c2ps[key])
end

"Finds the key of the cluster with the nearest center."
function nearest_other_center(c::ClusterSol, key)
    best_sqd = Inf
    best_key = 0
    for (key2, center) in c.centers
        if key2 != key
            sqd = sqdist(center, c.centers[key])
            if sqd < best_sqd
                best_sqd = sqd
                best_key = key2
            end
        end
    end
    return best_key
end

"Reassigns all points from the cluster with `key2` to the cluster with `key1`."
function merge_clusters!(c::ClusterSol, key1, key2)
    # Move assignment from key2 to key1 for all points in cluster with key2.
    for p in c.c2ps[key2]
        pop!(c.p2cs[p], key2)
        push!(c.p2cs[p], key1)
    end
    # Update the reverse index for the cluster with key1
    c.c2ps[key1] = union(c.c2ps[key1], c.c2ps[key2])
    # recompute the center for the cluster with key1
    update_center!(c, key1)
    # remove the cluster with key2
    delete!(c.c2ps, key2)
    delete!(c.centers, key2)
end

"Divides a cluster into `num_splits` equally sized clusters"
function split_cluster!(c::ClusterSol, key, num_splits)
    # Generate num_splits new keys
    new_keys = [rand(Int) for i in 1:num_splits]
    # Create the reverse indices for the new keys
    for new_key in new_keys
        c.c2ps[new_key] = Set{Int}()
    end
    # Reassign each point to one of the new clusters
    for (i,p) in enumerate(c.c2ps[key])
        new_key = new_keys[(i-1) % num_splits + 1]
        pop!(c.p2cs[p], key)
        push!(c.p2cs[p], new_key)
        push!(c.c2ps[new_key], p)
    end
    # Compute the centers for each new cluser
    for new_key in new_keys
        update_center!(c, new_key)
    end
    # Remove the old cluster
    delete!(c.c2ps, key)
    delete!(c.centers, key)
end

"""
Uses merging and splitting heuristics to ensure that each cluster satisfies the
capacity constraints. The resulting clustering might not have exactly `k`
clusters.
"""
function balancing_heuristics!(c::ClusterSol, l, L)
    # Merge small clusters with nearest cluster untill all clusters exceed l
    # points
    n = num_points(c)
    ks= cluster_keys(c)
    for key in ks
        if haskey(c.centers, key)
            if length(c.c2ps[key]) < l
                nearest_key = nearest_other_center(c, key)
                merge_clusters!(c, nearest_key, key)
            end
        end
    end
    # Split large clusters so that no cluster exceeds L points
    ks = cluster_keys(c)
    for key in ks
        if length(c.c2ps[key]) > L
            num_splits = ceil(Int, length(c.c2ps[key]) / L)
            split_cluster!(c, key, num_splits)
        end
    end
    for key in cluster_keys(c)
        @assert l <= length(c.c2ps[key]) <= L
    end
end

"""
Runs the kp-means algorithm and applies heuristics to enforce capacity
constraints on the cluster sizes
"""
function kpmeans(X, k, p, l, L)
    if p == 1 # If p = 1, use the nicely optimized Clustering.kmeans
        if k == 1
            centers = Dict{Int, Vector{Float64}}()
            centers[1] = squeeze(mean(X,2),2)
            p2cs = [Set{Int}(1) for i in 1:size(X,2)]
            c2ps = Dict{Int, Set{Int}}()
            c2ps[1] = Set{Int}(1:size(X,2))
            c = ClusterSol(X, centers, p2cs, c2ps)
        else
            c = ClusterSol(Clustering.kmeans(X, k), X)
        end
    else
        c = kpmeans(X, k, p)
    end
    balancing_heuristics!(c, l, L)
    return c
end
