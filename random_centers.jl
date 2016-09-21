#make each point its own cluster

function random_centers(X, k, p, l, L)
    @assert p == 1
    @assert k == size(X, 2)
    centers = Dict{Int, Vector{Float64}}()
    p2cs = [Set{Int}(i) for i in 1:size(X, 2)]
    c2ps = Dict{Int, Set{Int}}()
    for i in 1:size(X, 2)
        centers[i] = X[:, i]
        c2ps[i] = Set{Int}(i)
    end
    c = ClusterSol(X, centers, p2cs, c2ps)
    return c
end
