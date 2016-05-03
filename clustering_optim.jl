using JuMP, Gurobi

function solve_clustering_optimization(metric, k, p, ℓ, L;
                                       obj=:KMedian, solve_ip=false,
                                       threads=1, verbose=false)
    @assert obj in [:KMedian, :KMeans]
    N = size(metric)
    model = Model(solver = GurobiSolver(LogToConsole = verbose?1:0,
                                        Threads = threads))
    # Optimization variables
    if solve_ip
        @defVar(model, x[1:N, 1:N], Bin)
        @defVar(model, y[1:N], Bin)
    else
        @defVar(model, 0 <= x[1:N, 1:N] <= 1)
        @defVar(model, 0 <= y[1:N] <= 1)
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
    # Set the objective
    obj_expr = zero(AffExpr)
    for i in 1:N, j in 1:N
        if obj == :KMedian
            c = dist(metric, i, j)
        elseif obj == :KMeans
            c = dist(metric, i, j)^2
        end
        append!(obj_expr, c*x[i,j])
    end
    @setObjective(model, Min, obj_expr)
    # Sovle the model
    status = solve(model)

    return getValue(x), getValue(y)
end
