"Draw a sample from a discrete distribution."
function discrete_sample(ws, total)
    p = rand() * total
    t = 0.0
    for i in eachindex(ws)
        t += ws[i]
        if t > p
            return i
        end
    end
    error("discrete_sample failed")
    return -1
end

"Computes the squared distance between colum `aj` of `A` and colum `bj` of `B`."
function col_sqdist(A, aj, B, bj)
    sqdist = 0.0
    for i in 1:size(A,1)
        δ = A[i,aj] - B[i, bj]
        sqdist += δ*δ
    end
    sqdist
end

"Adds column `bj` of B into column `aj` of A."
function add_cols!(A, aj, B, bj)
    for i in 1:size(A,1)
        A[i,aj] += B[i,bj]
    end
end

"Scales column `j` of `A` by `s`"
function scale_col!(A, j, s)
    for i in 1:size(A,1)
        A[i,j] *= s
    end
end

"Copies column `bj` of `B` into column `aj` of `A`."
function copy_col!(A, aj, B, bj)
    for i in 1:size(A,1)
        A[i,aj] = B[i, bj]
    end
end

"Compute the squared distance between `u` and `v`."
function sqdist(u, v)
    total = 0.0
    for i in eachindex(u)
        δ = u[i] - v[i]
        total += δ*δ
    end
    return total
end
