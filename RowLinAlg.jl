module RowLinAlg

export sqnorm, sqdist, dist

using RowBlocks
import Base.dot, Base.norm

function dot(v::Row, u::Row)
    v_idxs = row_indices(v)
    v_vals = row_values(v)
    u_idxs = row_indices(u)
    u_vals = row_values(u)
    vi = ui = 1
    total = 0.0
    while vi <= length(v_idxs) && ui <= length(u_idxs)
        if v_idxs[vi] == u_idxs[ui]
            total += v_vals[vi] * u_vals[ui]
            vi += 1
            ui += 1
        elseif v_idxs[vi] < u_idxs[ui]
            vi += 1
        else
            ui += 1
        end
    end
    return total
end

sqnorm(v::Row) = dot(v,v)
norm(v::Row) = sqrt(sqnorm(v))
sqdist(v::Row, u::Row) = dot(v,v) + dot(u,u) - 2*dot(v,u)
dist(v::Row, u::Row) = sqrt(sqdist(v,u))

end
