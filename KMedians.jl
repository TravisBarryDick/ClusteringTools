module KMedians

export SparseLPSolution, exact_kmedian, kmedian_rounded_y, objective_value
export make_assignment_matrix

using Metrics
import Base.size

#################################################
# Code for solving the k-medians linear program #
#################################################

using JuMP, Clp, Cbc

function solve_kmedian_lp(metric, k, p, ℓ, L;
                          last_facility = -1, first_location = -1,
                          solve_ip = false, verbose = false, soft_capacities = false)

    N = size(metric)
    last_facility = last_facility == -1 ? N : last_facility
    num_facilities = last_facility

    first_location = first_location == -1 ? 1 : first_location
    num_locations = N - first_location + 1

    # Conversion from facility indices to metric indices (this is trivial for now)
    f2m(f) = f
    m2f(m) = m
    # Conversion from location indices to metric indices
    l2m(l) = l + first_location - 1
    m2l(m) = m - first_location + 1

    #model = Model(solver = GurobiSolver(LogToConsole = verbose ? 1 : 0, Threads=4))
    if solve_ip
        model = Model(solver = CbcSolver(LogLevel = verbose ? 4 : 0))
    else
        model = Model(solver = ClpSolver(LogLevel = verbose ? 4 : 0))
    end

    if solve_ip
        if soft_capacities
            @defVar(model, x[1:num_facilities, 1:num_locations] >= 0, Int)
            @defVar(model, y[1:num_facilities] >= 0, Int)
        else
            @defVar(model, x[1:num_facilities, 1:num_locations], Bin)
            @defVar(model, y[1:num_facilities], Bin)
        end
    else
        if soft_capacities
            @defVar(model, x[1:num_facilities, 1:num_locations] >= 0)
            @defVar(model, y[1:num_facilities] >= 0)
        else
            @defVar(model, 0 <= x[1:num_facilities, 1:num_locations] <= 1)
            @defVar(model, 0 <= y[1:num_facilities] <= 1)
        end
    end
    # Replication constraints
    for j in 1:num_locations
        @addConstraint(model, sum{x[i,j], i=1:num_facilities} == p)
    end
    # Capacity constraints
    for i in 1:last_facility
        @addConstraint(model, sum{x[i,j], j=1:num_locations} >= ℓ*y[i])
        @addConstraint(model, sum{x[i,j], j=1:num_locations} <= L*y[i])
    end
    # Number of centers
    @addConstraint(model, sum{y[i], i=1:num_facilities} == k)
    # Assignment constraints
    for i in 1:last_facility, j in 1:num_locations
        @addConstraint(model, x[i,j] >= 0)
        @addConstraint(model, x[i,j] <= y[i])
    end
    @setObjective(model, Min, sum{dist(metric,f2m(i),l2m(j))*x[i,j], i=1:num_facilities, j=1:num_locations})
    status = solve(model)
    return SparseLPSolution(getValue(x),
                            getValue(y),
                            N, last_facility, first_location)
end

####################
# Sparse Solutions #
####################

type SparseLPSolution
    centers::Dict{Int, Float64}
    assignments::Vector{Dict{Int, Float64}}
    clusters::Vector{Set{Int}}
    to_be_deleted::Set{Tuple{Int,Int}}
    eps::Float64
end

SparseLPSolution(N) = SparseLPSolution(Dict{Int,Float64}(),
                                       [Dict{Int,Float64}() for i in 1:N],
                                       [Set{Int}() for i in 1:N],
                                       Set{Tuple{Int,Int}}(), 1e-20)

function SparseLPSolution(x, y, N, last_facility, first_location)
    f2m(f) = f
    l2m(l) = l + first_location - 1
    sol = SparseLPSolution(N)
    for i in 1:last_facility
        set_opening!(sol, f2m(i), y[f2m(i)])
        for j in 1:(N-first_location+1)
            set_assignment!(sol, l2m(j), f2m(i), x[i,j])
        end
    end
    return sol
end

SparseLPSolution(x, y, N) = SparseLPSolution(x,y,N,N,1)

function make_assignment_matrix(s::SparseLPSolution, n, k, p)
    cs = Array(Int, k)
    as = Array(Int, p, n)

    cix = 1
    a_count = ones(Int, n)
    for (c,w) in s.centers
        copies = round(Int, w)
        for m in 1:copies
            cs[cix + m - 1] = c
        end

        copy_index = cix
        for x in s.clusters[c]
            for m in 1:round(Int, s.assignments[x][c])
                as[a_count[x], x] = copy_index
                a_count[x] += 1
                copy_index += 1
                if copy_index >= copies + cix
                    copy_index = cix
                end
            end
        end

        cix += copies
    end

    return as, cs
end

function objective_value(m::FiniteMetric, s::SparseLPSolution)
    cost = 0.0
    for i in 1:size(m)
        for (c, w) in s.assignments[i]
            cost += dist(m, i, c) * w
        end
    end
    return cost
end

function objective_value(data, s::SparseLPSolution)
    return objective_value(precompute_metric(data), s)
end

function mean_distance(m::FiniteMetric, s::SparseLPSolution, i)
    p = sum(values(s.assignments[i]))
    cost = 0.0
    for (c, w) in s.assignments[i]
        cost += w*dist(m, i, c)
    end
    return cost / p
end

function set_assignment!(s::SparseLPSolution, i, j, w)
    if w ≥ s.eps
        s.assignments[i][j] = w
        push!(s.clusters[j], i)
    else
        delete!(s.assignments[i], j)
        if i in s.clusters[j]
            push!(s.to_be_deleted, (i,j))
        end
    end
end

function clean_clusters!(s::SparseLPSolution)
    for (i,j) in s.to_be_deleted
        delete!(s.clusters[j], i)
    end
    empty!(s.to_be_deleted)
end

function set_opening!(s::SparseLPSolution, j, w)
    if w ≥ s.eps
        s.centers[j] = w
    else
        delete!(s.centers, j)
    end
end

# Moves δ opening from center at i to center at j
function move_y!(s::SparseLPSolution, i, j, δ, eps=1e-8)
    γ = δ / get(s.centers, i, 0.0)
    set_opening!(s, i, get(s.centers, i, 0.0) - δ)
    set_opening!(s, j, get(s.centers, j, 0.0) + δ)
    for p in s.clusters[i]
        Δx = get(s.assignments[p], i, 0.0) * γ
        set_assignment!(s, p, i, get(s.assignments[p], i, 0.0) - Δx)
        set_assignment!(s, p, j, get(s.assignments[p], j, 0.0) + Δx)
    end
    clean_clusters!(s)
end

######################
# Rounding Procedure #
######################

function monarch_dist(m::FiniteMetric, monarchs, i)
    min_dist = Inf # best distance so far
    best_mix = 0   # corresponding monarch index
    for mix in 1:length(monarchs)
        d = dist(m, monarchs[mix], i)
        if d < min_dist
            min_dist = d
            best_mix = mix
        end
    end
    return min_dist, best_mix
end

### Monarch Procedure

function monarch_procedure(m::FiniteMetric, s::SparseLPSolution)
    N = size(m)
    # Compute the set of monarchs
    C = [mean_distance(m, s, i) for i in 1:N]
    perm = sortperm(C)
    monarchs = Array(Int, 0)
    for ix in perm
        d, mix = monarch_dist(m, monarchs, ix)
        if d > 4*C[ix]
            push!(monarchs, ix)
        end
    end
    # Compute the empire for each monarch
    empires = [Array(Int,0) for monarch in monarchs]
    for ix in 1:N
        d, mix = monarch_dist(m, monarchs, ix)
        push!(empires[mix], ix)
    end
    return monarchs, empires
end

function check_monarch_properties(m::FiniteMetric, s::SparseLPSolution, p,
                                  monarchs, empires)
    N = size(m)
    N′ = length(monarchs)
    C = [mean_distance(m, s, i) for i in 1:N]
    # ∀ l ∈ N′, ∀ j ∈ M_l, dist(m, l, j) < 4*C[j]
    for (l, e) in zip(monarchs, empires)
        for j in e
            @assert dist(m, j, l) <= 4*C[j]
        end
    end
    # ∀ l ≠ l′ ∈ N′, dist(m, l, l′) > 4*max(C[l], C[l′])
    # ∀ l ≠ l′ ∈ N′, M_l ∩ M_l′ = ∅
    for (l,e) in zip(monarchs, empires), (l′,e′) in zip(monarchs, empires)
        if l == l′; continue; end
        @assert dist(m, l, l′) > 4*max(C[l], C[l′])
        @assert length(e ∩ e′) == 0
    end
    # ∪_{l ∈ N′} M_l = N
    seen = zeros(Bool, N)
    for e in empires, i in e
        seen[i] = true
    end
    @assert all(seen)
    # ∀ l ∈ N′, ∑_{j ∈ M_l} y_j ≥ p/2
    for (l, e) in zip(monarchs, empires)
        total_opening = 0.0
        for ix in e
            total_opening += get(s.centers, ix, 0.0)
        end
        @show total_opening
        @assert total_opening >= p/2
    end
end

function round_ys!(m::FiniteMetric, s::SparseLPSolution, monarchs, empires)
    for (monarch, empire) in zip(monarchs, empires)
        # Calculate the total opening in this empire
        Y = 0.0
        for i in empire
            Y += get(s.centers, i, 0.0)
        end
        # Calculate the target opening
        z = Y / floor(Y)
        # Sort the empire by distance to the monarch
        sort!(empire, by = i -> dist(m, monarch, i))
        inner_ix = 1
        for p in length(empire):-1:(floor(Int, Y)+1)
            source = empire[p]
            while get(s.centers, source, 0.0) > 0.0
                dest = empire[inner_ix]
                to_move = get(s.centers, source, 0.0)
                space = z - get(s.centers, dest, 0.0)
                amount = min(to_move, space)
                move_y!(s, source, dest, amount)
                if get(s.centers, dest, 0.0) >= z
                    inner_ix += 1
                end
            end
        end
    end
    for (c, w) in s.centers
        set_opening!(s, c, 1.0)
    end
end

#####################
# Wrapper Functions #
#####################

function exact_kmedian(data, k, p, ℓ, L; kwargs...)
    metric = precompute_metric(data)
    return exact_kmedian(metric, k, p, ℓ, L; kwargs...)
end

function exact_kmedian(data, k; kwargs...)
    exact_kmedian(data, k, 1, 0, size(data,2); kwargs...)
end

function exact_kmedian(metric::FiniteMetric, k, p, ℓ, L; kwargs...)
    solve_kmedian_lp(metric, k, p, ℓ, L, solve_ip=true; kwargs...)
end

function exact_kmedian(metric::FiniteMetric, k; kwargs...)
    exact_kmedian(metric, k, 1, 0, size(metric); kwargs...)
end

function kmedian_rounded_y(metric::FiniteMetric, k, p, ℓ, L; kwargs...)
    @assert p ≥ 2
    sol = solve_kmedian_lp(metric, k, p, ℓ, L; kwargs...)
    monarchs, empires = monarch_procedure(metric, sol)
    round_ys!(metric, sol, monarchs, empires)
    return sol
end

function kmedian_rounded_y(data, k, p, ℓ, L; verbose = false)
    metric = precompute_metric(data)
    kmedian_rounded_y(metric, k, p, ℓ, L; verbose = verbose)
end


include("round_xs.jl")

end
