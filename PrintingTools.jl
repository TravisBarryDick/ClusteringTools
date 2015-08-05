module PrintingTools

export print_clustering

function print_clustering(as, cs)
    p, n = size(as)
    k = length(cs)
    clusters = [Set{Int}() for i in 1:k]
    for j in 1:n
        for i in 1:p
            push!(clusters[as[i,j]], j)
        end
    end

    colors = [:red, :yellow, :magenta, :green, :cyan, :blue]
    for i in 1:size(as,2)
        print('.')
    end
    println()
    for i in 1:k
        color = colors[(i-1)%length(colors) + 1]
        for p in 1:n
            if p == cs[i]
                if cs[i] in clusters[i]
                    print_with_color(color, "x")
                else
                    print("o")
                end
            elseif p in clusters[i]
                print_with_color(color, "-")
            else
                print(' ')
            end
        end
        println()
    end
end

end
