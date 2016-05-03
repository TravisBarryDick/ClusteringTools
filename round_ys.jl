"""
Moves δ opening from center `i` to center `j` (and proportionally moves
assignments).
"""
function move_y!(s::SparseLPSolution, i, j, δ)
    γ = δ / get_y(s, i)
    set_y!(s, i, get_y(s, i) - δ)
    set_y!(s, j, get_y(s, j) + δ)
    for p in get_cluster(s, i)
        Δx = get_x(s, i, p) * γ
        set_x!(s, i, p, get_x(s, i, p) - Δx)
        set_x!(s, j, p, get_x(s, j, p) + Δx)
    end
    sparsify!(s)
end

"Computes the mean cost of point `j` to its assigned centers."
function mean_cost(m::FiniteMetric, s::SparseLPSolution, j, obj)
    cost = 0.0
    p = 0.0
    for (i, xij) in s.x[j]
        p += xij
        if obj == :KMedian
            cost += dist(m, i, j) * xij
        elseif obj == :KMeans
            cost += dist(m, i, j)^2 * xij
        end
    end
    return cost / p
end

function monarch_procedure(m::FiniteMetric, s::SparseLPSolution, obj)
    if obj == :KMedian
        ρ = 2
    elseif obj == :KMeans
        ρ = 4
    end

    N = size(m)
    C = [mean_cost(m, s, j, obj) for j in 1:N]
    perm = sortperm(C)
    monarchs = Array(Int, 0)

    function mdist(ix)
        min_dist = Inf
        best_mix = 0
        for mix in 1:length(monarchs)
            d = dist(m, monarchs[mix], ix)
            if d < min_dist
                min_dist = d
                best_mix = mix
            end
        end
        return min_dist, best_mix
    end

    # Compute the set of monarchs
    for ix in perm
        d, mix = mdist(ix)
        if d > 2*ρ*C[ix]
            push!(monarchs, ix)
        end
    end

    # Compute the empire for each monarch
    empires = [Array(Int, 0) for monarch in monarchs]
    for ix in 1:N
        d, mix = mdist(ix)
        push!(empires[mix], ix)
    end

    return monarchs, empires
end

function round_ys!(m::FiniteMetric, s::SparseLPSolution, obj)
    monarchs, empires = monarch_procedure(m, s, obj)
    for (monarch, empire) in zip(monarchs, empires)
        Y = 0.0
        for i in empire
            Y += get_y(s, i)
        end
        z = Y / floor(Y)
        sort!(empire, by = i -> dist(m, monarch, i))
        inner_ix = 1
        for p in length(empire):-1:(floor(Int, Y)+1)
            source = empire[p]
            while get_y(s, source) > sparse_eps
                dest = empire[inner_ix]
                to_move = get_y(s, source)
                space = z - get_y(s, dest)
                amount = min(to_move, space)
                move_y!(s, source, dest, amount)
                if get_y(s, dest) >= z
                    inner_ix += 1
                end
            end
            set_y!(s, source, 0.0)
            sparsify!(s)
        end
    end
    for c in keys(s.y)
        set_y!(s, c, 1.0)
    end
    sparsify!(s)
    nothing
end
