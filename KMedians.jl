module KMedians

export SparseLPSolution, exact_kmedian

using Metrics
import Base.size

#################################################
# Code for solving the k-medians linear program #
#################################################

using JuMP, Gurobi

function solve_kmedian_lp(metric, k, p, ℓ, L; solve_ip=false, verbose=false)
    N = size(metric)
    @assert ℓ*k <= p*N "Lower capacity can't be satisfied."
    @assert L*k >= p*N "Upper capacity can't be satisfied."

    # Use JuMP to build the LP
    model = Model(solver = GurobiSolver(LogToConsole= verbose ? 1 : 0, Threads=4))
    # Defining the variables
    if solve_ip
        @defVar(model, x[1:N, 1:N], Bin)
        @defVar(model, y[1:N], Bin)
    else
        @defVar(model, 0 <= x[1:N, 1:N] <= 1) # assignment variables
        @defVar(model, 0 <= y[1:N] <= 1) # opening variables
    end
    # Replication constraints
    for j in 1:N
        @addConstraint(model, sum{x[i,j], i=1:N} == p)
    end
    # Capacity constraints
    for i in 1:N
        @addConstraint(model, sum{x[i,j], j=1:N} >= ℓ*y[i])
        @addConstraint(model, sum{x[i,j], j=1:N} <= L*y[i])
    end
    # Number of centers
    @addConstraint(model, sum{y[i], i=1:N} == k)
    # Assignment constraints
    for i in 1:N, j in 1:N
        @addConstraint(model, x[i,j] >= 0)
        @addConstraint(model, x[i,j] <= y[i])
    end
    @setObjective(model, Min, sum{dist(metric,i,j)*x[i,j], i=1:N, j=1:N})
    status = solve(model)
    return SparseLPSolution(getValue(x).innerArray,
                            getValue(y).innerArray)
end

####################
# Sparse Solutions #
####################

type SparseLPSolution
    centers::Dict{Int, Float64}
    assignments::Vector{Dict{Int, Float64}}
    clusters::Dict{Int, Set{Int}}
    eps::Float64
end
#=
Invariants for SparseLPSolution instances:
1. No zero centers:
    ∀ j : haskey(centers, j) ⇒ centers[j] > eps
2. No zero assignments:
    ∀ i, ∀ j : haskey(assignments[i], j) ⇒ assignments[i][j] > eps
3. Clusters and centers are in one-to-one correspondence always
    ∀ i : haskey(centers, j) ⇔ haskey(clusters, j)
4. Can't over-assign
    ∀ i, ∀ j : assignments[i][j] ≤ centers[j]
5. Assignments must be consistent with clusters
    ∀ i, ∀ j : haskey(assignments[i], j) ⇒ i ∈ clusters[j]
=#

SparseLPSolution(N) = SparseLPSolution(Dict{Int,Float64}(),
                                       [Dict{Int,Float64}() for i in 1:N],
                                       Dict{Int, Set{Int}}(), 1e-10)

function objective_value(m::FiniteMetric, s::SparseLPSolution)
    cost = 0.0
    for i in 1:size(m)
        for (c, w) in s.assignments[i]
            cost += dist(m, i, c) * w
        end
    end
    return cost
end

function mean_distance(m::FiniteMetric, s::SparseLPSolution, i)
    cost = 0.0
    for (c, w) in s.assignments[i]
        cost += w*dist(m, i, c)
    end
    return cost / length(s.assignments[i])
end

function set_assignment!(s::SparseLPSolution, i, j, w)
    if abs(w) < s.eps
        delete!(s.assignments[i], j)
        if haskey(s.clusters, j)
            delete!(s.clusters[j], i)
        end
    else
        @assert haskey(s.centers, j) && haskey(s.clusters, j) && s.centers[j] ≥ w
        s.assignments[i][j] = w
        push!(s.clusters[j], i)
    end
end

function set_opening!(s::SparseLPSolution, j, w)
    if abs(w) < s.eps
        @assert !haskey(s.clusters, j) || length(s.clusters[j]) == 0
        delete!(s.centers, j)
        delete!(s.clusters, j)
    else
        s.centers[j] = w
        s.clusters[j] = Set{Int}()
    end
end

function SparseLPSolution(x, y)
    N = size(x,1)
    s = SparseLPSolution(N)
    for j in 1:N
        set_opening!(s, j, y[j])
    end
    for i in 1:N, j in 1:N
        set_assignment!(s, i, j, x[j,i])
    end
    return s
end

# Moves δ opening from center at i to center at j
function move_y!(s::SparseLPSolution, i, j, δ, eps=1e-8)
    γ = δ / get(s.centers, i, 0.0) # fraction of x assignments to move
    set_opening!(s, i, get(s.centers, i, 0.0) - δ)
    set_opening!(s, j, get(s.centers, j, 0.0) + δ)
    if haskey(s.centers, i)
        for p in s.clusters[i]
            Δx = get(s.assignments[p], i, 0.0) * γ
            set_assignment!(s, p, i, get(s.assignments[p], i, 0.0) - Δx)
            set_assignment!(s, p, j, get(s.assignments[p], i, 0.0) + Δx)
        end
    end
end

#####################
# Wrapper Functions #
#####################

function exact_kmedian(metric, k, p, ℓ, L; verbose=false)
    solve_kmedian_lp(metric, k, p, ℓ, L, solve_ip=true, verbose=verbose)
end

function exact_kmedian(metric, k; verbose=false)
    exact_kmedian(metric, k, 1, 0, size(metric); verbose=verbose)
end

end
