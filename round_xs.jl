export kmedian_round

function output_dimacs(filename, metric::FiniteMetric, sol::SparseLPSolution, p, ℓ, L)
    centers = collect(keys(sol.y))
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
    centers = collect(keys(sol.y))
    k = length(centers)
    new_sol = SparseLPSolution(N)
    for c in centers
        set_y!(new_sol, c, 1.0)
    end
    f = open(filename)
    for line in eachline(f)
        flow = map( x -> parse(Int, strip(x)), split(line))
        if (length(flow) != 3) ; continue; end
        if (flow[1] > N) ; continue ; end #ignore edges to sink
        @assert ( flow[2] > N && flow[2] <= N+k )
        set_x!(new_sol, centers[flow[2] - N], flow[1], flow[3])
    end
    return new_sol
end

function round_xs(metric::FiniteMetric, sol::SparseLPSolution, p, ℓ, L;
    verbose = false)
    infile = "./$(hash(rand()))"
    outfile = "./$(hash(rand()))"
    output_dimacs(infile, metric, sol, p, ℓ, L)
    # to compile: g++ solve_mcf.cc -I ../lemon/include/ -L ../lemon/lib/ -lemon -o mcf
    mcf_path = joinpath(dirname(@__FILE__()), "mcf")
    run(`$mcf_path $infile $outfile`)
    N = size(metric)
    new_sol = parse_mcf_output(outfile, sol,  N)
    rm(infile)
    rm(outfile)
    sparsify!(new_sol)
    return new_sol
end
