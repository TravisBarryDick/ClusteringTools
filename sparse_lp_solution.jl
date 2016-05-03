"""
An easy-to-modify representation of a sparse solution to the clustering LP. Note
that the x variable gives fast indexing from points to the centers they are
assigned to, while c gives fast indexing from centers to the points assigned to
them.
"""
type SparseLPSolution
    y::Dict{Int, Float64} # get(y, i, 0) = opening of point i
    x::Vector{Dict{Int, Float64}} # get(x[j], i, 0) = assignment of j to i
    c::Vector{Dict{Int, Float64}} # get(c[i], j, 0) = assignment of j to i
    # Some extra state to keep things sparse
    _xz::Set{Tuple{Int,Int}} # assignment variables equal to zero
    _yz::Set{Int} # opening variables equal to zero
end

const sparse_eps = 1e-8

"Constructs the all-zero sparse LP solution."
function SparseLPSolution(N::Int)
    y = Dict{Int, Float64}()
    x = [Dict{Int, Float64}() for i in 1:N]
    c = [Dict{Int, Float64}() for i in 1:N]
    _xz = Set{Tuple{Int,Int}}()
    _yz = Set{Int}()
    return SparseLPSolution(y, x, c, _xz, _yz)
end

"Constructs a `SparseLPSolution` from dense matrices x and y."
function SparseLPSolution(x::Array{Float64,2}, y::Vector{Float64})
    N = size(y, 1)
    sol = SparseLPSolution(N)
    for i in 1:N
        set_y!(sol, i, y[i])
        for j in 1:N
            set_x!(sol, i, j, x[i,j])
        end
    end
    return sol
end

"Set the assignment of node `j` to node `i` to be `v`."
function set_x!(s::SparseLPSolution, i, j, v)
    if abs(v) >= sparse_eps || haskey(s.x[j], i)
        s.x[j][i] = v
        s.c[i][j] = v
        if abs(v) < sparse_eps
            push!(s._xz, (i,j))
        end
    end
    nothing
end

"Returns the assignment of node `j` to node `i`."
get_x(s::SparseLPSolution, i, j) = get(s.x[j], i, 0.0)

"Set the opening of node `i` to `v`."
function set_y!(s::SparseLPSolution, i, v)
    if abs(v) >= sparse_eps || haskey(s.y, i)
        s.y[i] = v
        if abs(v) < sparse_eps
            push!(s._yz, i)
        end
    end
    nothing
end

"Returns the opening of node `i`."
get_y(s::SparseLPSolution, i) = get(s.y, i, 0.0)

"Returns the indices of points with possibly non-zero assignment to center `i`."
get_cluster(s::SparseLPSolution, i) = keys(s.c[i])

"Returns the indices of centers that `j` has possibly non-zero assignment to."
get_centers(s::SparseLPSolution, j) = keys(s.x[j])

"Returns the indices of all centers in the clustering."
get_centers(s::SparseLPSolution) = keys(s.y)

"Deletes all variables whose value is < sparse_eps."
function sparsify!(s::SparseLPSolution)
    for (i,j) in s._xz
        if abs(s.x[j][i]) < sparse_eps
            delete!(s.x[j], i)
        end
        if abs(s.c[i][j]) < sparse_eps
            delete!(s.c[i], j)
        end
    end
    empty!(s._xz)
    for i in s._yz
        if abs(s.y[i]) < sparse_eps
            delete!(s.y, i)
        end
    end
    empty!(s._yz)
end

function objective_value(m::FiniteMetric, s::SparseLPSolution, obj)
    total = 0.0
    for i in get_centers(s)
        for (j,xij) in s.c[i]
            if obj == :KMedian
                total += dist(m, i, j) * xij
            elseif obj == :KMeans
                total += dist(m, i, j)^2 * xij
            end
        end
    end
    return total
end

function ClusterSol(s::SparseLPSolution, X::Array{Float64, 2})
    n = size(X,2)
    centers = [k => X[:,k] for k in keys(s.centers)]
    p2cs = [Set(collect(keys(as))) for as in s.assignments]
    c2ps = [k => s.clusters[k] for k in keys(s.centers)]
    return ClusterSolution(X, centers, p2cs, c2ps)
end
