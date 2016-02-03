function round_xs(metric, centers, p, ℓ, L;
                          verbose = false)
#centers is an array with the indices of centers

    N = size(metric)
    k = length(centers)
    new_L = ceil(1.0 * L * (p+2) / p)

    model = Model(solver = GurobiSolver(LogToConsole = verbose ? 1 : 0, Threads=4))
  
    # Define Variables
    @defVar(model, x[1:k,  1:N], Bin)
    # Replication constraints
    for j in 1:N
        @addConstraint(model, p/2 <= sum{x[i,j], i=1:k} <= p)
    end
    # Capacity constraints
    for i in 1:k
        @addConstraint(model, sum{x[i,j], j=1:N} >= ℓ)
        @addConstraint(model, sum{x[i,j], j=1:N} <= new_L)
    end
    @setObjective(model, :Min, sum{dist(metric,centers(i),j)*x[i,j], i=1:k, j=1:N})
    status = solve(model)
    return getValue(x)
end
