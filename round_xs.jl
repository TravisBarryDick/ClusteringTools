export kmedian_round

function round_xs(metric::FiniteMetric, sol::SparseLPSolution, p, ℓ, L;
        verbose = false)
        #centers is an array with the indices of centers
        centers = collect(keys(sol.centers))

        N = size(metric)
        k = length(centers)
        new_L = ceil(1.0 * L * (p+2) / p)

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
			println("Assign center $(centers[i]) to point $(j)")
		end
	end
	return new_sol
end

function kmedian_round(metric::FiniteMetric, k, p, ℓ, L; kwargs...)
        @assert p ≥ 2 
        sol = solve_kmedian_lp(metric, k, p, ℓ, L; kwargs...)
        monarchs, empires = monarch_procedure(metric, sol)
        round_ys!(metric, sol, monarchs, empires)
        new_sol = round_xs(metric, sol, p, ℓ, L, verbose = true)
        return new_sol
end

function kmedian_round(data, k, p, ℓ, L; kwargs...)
        metric = precompute_metric(data)
        return kmedian_round(metric, k, p, ℓ, L; kwargs...)
end
