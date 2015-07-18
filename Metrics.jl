module Metrics

export FiniteMetric, dist
export PrecomputedMetric, precompute_metric
export ImplicitMetric, make_implicit_metric
import Base.size

using RowBlocks, RowLinAlg

abstract FiniteMetric
size(m::FiniteMetric) = error("Must implement size for your metric!")
dist(m::FiniteMetric, i, j) = error("Must implement dist for your metric!")

type PrecomputedMetric <: FiniteMetric
    d::Array{Float64, 2}
end
size(m::PrecomputedMetric) = size(m.d,1)
dist(m::PrecomputedMetric, i, j) = m.d[i,j]

function precompute_metric(xs::Array{Float64, 2})
    n = size(xs, 2)
    d = Array(Float64, n, n)
    for i in 1:n, j in 1:i
        distance = dot(xs[:,i], xs[:,i]) + dot(xs[:,j], xs[:,j]) -
                   2*dot(xs[:,i], xs[:,j])
        d[i,j] = distance
        d[j,i] = distance
    end
    return PrecomputedMetric(d)
end

function precompute_metric(xs)
    n = length(xs)
    d = Array(Float64, n, n)
    for i in 1:n, j in 1:i
        distance = dot(xs[i],xs[i]) + dot(xs[j],xs[j]) - 2*dot(xs[i],xs[j])
        d[i,j] = distance
        d[j,i] = distance
    end
    return PrecomputedMetric(d)
end

end
