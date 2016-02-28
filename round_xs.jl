export kmedian_round

function round_xs_ip(metric::FiniteMetric, sol::SparseLPSolution, p, ℓ, L;
    verbose = false)
    #output_dimacs("./temp_out", metric, sol, p, ℓ, L)
    #centers is an array with the indices of centers
    centers = collect(keys(sol.centers))

    N = size(metric)
    k = length(centers)
    new_L = ceil(Int, 1.0 * L * (p+2) / p)

    model = Model(solver = GurobiSolver(LogToConsole = verbose ? 1 : 0, Threads=4))

    # Define Variables
    @defVar(model, x[1:k,  1:N], Bin)
    #@defVar(model, 0 <= x[1:k,  1:N] <= 1)
    # Replication constraints
    for j in 1:N
        @addConstraint(model, p/2 <= sum{x[i,j], i=1:k} <= p)
    end
    # Capacity constraints
    for i in 1:k
        @addConstraint(model, sum{x[i,j], j=1:N} >= ℓ)
        @addConstraint(model, sum{x[i,j], j=1:N} <= new_L)
    end
    @setObjective(model, :Min, sum{dist(metric,centers[i],j)*x[i,j], i=1:k, j=1:N})
    status = solve(model)
    xs = getValue(x)
    # output
    new_sol = SparseLPSolution(N)
    for c in centers
        set_opening!(new_sol, c, 1.0)
    end
    for i in 1:k, j in 1:N
        if (xs[i, j] >= 1 - new_sol.eps)
            set_assignment!(new_sol, j, centers[i], xs[i, j])
            #println("Assign center $(centers[i]) to point $(j)")
        end
    end
    return new_sol
end

function output_dimacs(filename, metric::FiniteMetric, sol::SparseLPSolution, p, ℓ, L)
    centers = collect(keys(sol.centers))
    N = size(metric)
    k = length(centers)
    new_L = ceil(Int, 1.0 * L * (p+2) / p)
    outfile = open(filename, "w")
    println(outfile, "c DIMACS format of graph for min cost flow")
    #problem description
    println(outfile, "p min $(N+k+1) $(k*(N+1))")

    #node descriptions: points, centers and sink
    for i in 1:N
        println(outfile, "n $(i) $(p)")
    end
    for j in 1:k
        println(outfile, "n $(j+N) $(-ℓ)")
    end
    println(outfile, "n $(N+k+1) $(-N*p + k*ℓ)")

    #edge descriptions: points-to-centers and centers-to-sink
    for i in 1:N
        for j in 1:k
            println(outfile, "a $(i) $(j+N) 0 2 $(dist(metric,centers[j],i))")
        end
    end
    for j in 1:k
        println(outfile, "a $(j+N) $(N+k+1) 0 $(new_L - ℓ) 0.0")
    end
    close(outfile)
end

function parse_mcf_output(filename, sol::SparseLPSolution, N)
    #each line of output is "from to flow"
    centers = collect(keys(sol.centers))
    k = length(centers)
    new_sol = SparseLPSolution(N)
    for c in centers
        set_opening!(new_sol, c, 1.0)
    end
    f = open(filename)
    for line in eachline(f)
        flow = map( x -> parse(Int, strip(x)), split(line))
        if (length(flow) != 3) ; continue; end
        if (flow[1] > N) ; continue ; end #ignore edges to sink
        @assert ( flow[2] > N && flow[2] <= N+k )
        set_assignment!(new_sol, flow[1], centers[flow[2] - N], flow[3])
    end
    return new_sol
end

function round_xs_lemon(metric::FiniteMetric, sol::SparseLPSolution, p, ℓ, L;
    verbose = false)
    infile = "./temp1"
    outfile = "./temp2"
    output_dimacs(infile, metric, sol, p, ℓ, L)
    # to compile: g++ solve_mcf.cc -I ../lemon/include/ -L ../lemon/lib/ -lemon -o mcf
    mcf_path = joinpath(dirname(@__FILE__()), "mcf")
    run(`$mcf_path $infile $outfile`)
    N = size(metric)
    new_sol = parse_mcf_output(outfile, sol,  N)
    run(`rm $infile $outfile`)
    return new_sol
end

function kmedian_round(metric::FiniteMetric, k, p, ℓ, L; kwargs...)
    @assert p ≥ 2
    sol = solve_kmedian_lp(metric, k, p, ℓ, L; kwargs...)
    monarchs, empires = monarch_procedure(metric, sol)
    round_ys!(metric, sol, monarchs, empires)
    new_sol = round_xs_lemon(metric, sol, p, ℓ, L, verbose = true)
    return new_sol
end

function kmedian_round(data, k, p, ℓ, L; kwargs...)
    metric = precompute_metric(data)
    return kmedian_round(metric, k, p, ℓ, L; kwargs...)
end
