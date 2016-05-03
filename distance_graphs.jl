type DistanceGraph
    num_nodes::Int
    neighbors::Vector{Set{Int}}
    bin_size::Vector{Int}
    start_idx::Vector{Int}
end

function DistanceGraph(bin_sizes::Vector{Int})
    num_nodes = length(bin_sizes)
    start_idx = cumsum([1; bin_sizes])[1:end-1]
    return DistanceGraph(num_nodes, [Set{Int}() for i in 1:num_nodes], bin_sizes, start_idx)
end

function add_edge!(sg::DistanceGraph, a, b)
    push!(sg.neighbors[a], b)
    push!(sg.neighbors[b], a)
end

function node_group(sg::DistanceGraph, a)
    return sg.start_idx[a]:(sg.start_idx[a] + sg.bin_size[a] - 1)
end

total_size(sg::DistanceGraph) = sum(sg.bin_size)

function make_metric(sg::DistanceGraph)
    n_pts = total_size(sg)
    M = 2*ones(n_pts, n_pts) - 2*eye(n_pts)
    for a in 1:sg.num_nodes
        for b in sg.neighbors[a]
            for ai in node_group(sg, a), bi in node_group(sg, b)
                M[ai, bi] = 1.0
            end
        end
    end
    return PrecomputedMetric(M)
end

function node_assignment(sg::DistanceGraph, sol::SparseLPSolution, a)
    counts = Dict{Int, Int}()
    for i in node_group(sg, a)
        for (c,w) in sol.assignments[i]
            counts[c] = get(counts, c, 0) + 1
        end
    end
    return counts
end
