#2-D example. Mixture of gaussian
# True class label is determined by (global separator + random noise)

function generate_syn_1(ngaussian, nppg, seed = 0) #ppg = points per gaussian
    rng = MersenneTwister(seed)
    points = randn(rng, ngaussian * nppg, 2)
    labels = zeros(ngaussian * nppg)
    for i in 1:ngaussian
        # means of cluster i = (3*i,0) ; covariance = I
        m = 3*i
        range = (1 + (i-1)*nppg) : i*nppg
        points[range, 1] += m
        theta = rand(rng) * 90 #angle with x-axis: between 0 and 90 deg
        slope = tan(theta * pi / 180)
        # equation of separator: y = slope * (x - m)
        labels[range] = sign( points[range, 2] - slope * (points[range, 1] - m))
    end
    return points, ceil(Int, labels)
end

