module RowBlocks

export RowBlock, labels, Row, row_values, row_indices
export parse_libsvm, save_libsvm

import ArrayViews
import Base.length, Base.getindex
import Base.start, Base.next, Base.done

type RowBlock{I,L,R}
    offset::Vector{Int}
    value::Vector{R}
    index::Vector{I}
    label::Vector{L}
    dim::Int
end

function RowBlock(I::Type, L::Type, R::Type)
    RowBlock(Int[1], Array(R, 0), Array(I, 0), Array(L, 0), 0)
end

length(block::RowBlock) = length(block.offset) - 1
dimension(block::RowBlock) = block.dim
getindex(block::RowBlock, i) = Row(block, i)

labels(block::RowBlock) = block.label

immutable Row{I,L,R}
    block::RowBlock{I,L,R}
    index::Int
end

length(r::Row) = r.block.offset[r.index + 1] - r.block.offset[r.index]
dimension(r::Row) = dimension(r.block)
getindex(r::Row, i) = (r.block.index[r.index + i - 1], r.block.value[r.index + i - 1])

# Expose an array view for the indices and the values for this row
for (fn, arr) in [(:row_indices, :index), (:row_values, :value)]
    @eval function $fn(r::Row)
        first_index = r.block.offset[r.index]
        last_index = r.block.offset[r.index + 1] - 1
        return ArrayViews.view(r.block.$arr, first_index:last_index)
    end
end

# Make both RowBlocks and Rows iterable
for t in [:RowBlock, :Row]
    @eval start(x::$t) = 1
    @eval next(x::$t, state) = x[state], state+1
    @eval done(x::$t, state) = state == length(x) + 1
end

function parse_libsvm_line(block::RowBlock, line)
    space_ix = findfirst(line, ' ')
    y = parse(Int, line[1:space_ix-1])
    push!(block.label, y)
    for token in split(strip(line[space_ix+1:end]), ' ')
        colon_ix = findfirst(token, ':')
        ix = parse(Int, token[1:colon_ix-1])
        val = parse(Float64, token[colon_ix+1:end])
        push!(block.index, ix)
        push!(block.value, val)
        block.dim = max(block.dim, ix)
    end
    push!(block.offset, length(block.value)+1)
end

function parse_libsvm(stream::IO)
    block = RowBlock(Int, Int, Float64)
    for line in eachline(stream)
        parse_libsvm_line(block, line)
    end
    return block
end

function parse_libsvm(stream::IO, filter)
    block = RowBlock(Int, Int, Float64)
    for (i, line) in enumerate(eachline(stream))
        if filter(i, line)
            parse_libsvm_line(block, line)
        end
    end
    return block
end

function parse_libsvm(path, filter = nothing)
    f = open(path, "r")
    if filter != nothing
        block = parse_libsvm(f, filter)
    else
        block = parse_libsvm(f)
    end
    close(f)
    return block
end

function save_libsvm(stream::IO, block::RowBlock)
    for i in 1:length(block)
        print(stream, block.label[i], ' ')
        row = block[i]
        for j in 1:length(row)
            idx, val = row[j]
            print(stream, idx, ':', val, j == length(row) ? '\n' : ' ')
        end
    end
end

function save_libsvm(path, block::RowBlock)
    f = open(path, "w")
    save_libsvm(f, block)
    close(f)
end

end
