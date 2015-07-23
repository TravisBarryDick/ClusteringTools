module PlottingTools

using PyPlot

export plot_clustering, plot_data

function plot_clustering(data, sol)
    for cix in keys(sol.centers)
        ixs = collect(sol.clusters[cix])
        plot(data[1,ixs]', data[2,ixs]', ".k", ms=5)
        for ix in ixs
            plot([data[1,ix], data[1,cix]],
                 [data[2,ix], data[2,cix]],
                 "-k", alpha=sol.assignments[ix][cix]/2)
        end
    end

    for cix in keys(sol.centers)
        plot(data[1,cix], data[2,cix], ".", ms=10)
    end
end

function plot_data(data)
    plot(data[1,:]', data[2,:]', "k.")
end

end
