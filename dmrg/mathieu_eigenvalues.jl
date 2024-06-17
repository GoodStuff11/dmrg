using LaTeXStrings
using MathieuFunctions

function interleave(a::Vector{T}, b::Vector{T}) where T
    length(a) == length(b) || error("length mismatch")
    c = Vector{T}(undef, 2*length(a))
    i = 0
    for (x, y) in zip(a, b)
        c[i += 1] = x
        c[i += 1] = y
    end
    return c
end

function filter_checkerboard(x)
    xlen, ylen = size(x)
    output = zeros(xlen*ylen ÷ 2)
    k = 1
    for i=1:xlen
        for j=1:ylen
            if (i÷2 + (j÷2))%2 == 0
                output[k] = x[i,j]
                k += 1
            end
        end
    end 
    return output
end

f(q) = interleave(charA(q;k=0:10), charB(q;k=1:11))

p = plot(xlabel=L"g", ylabel=L"E", legend=nothing)
g_values = 0:0.01:30
data = reduce(vcat,transpose.([collect(Iterators.flatten(filter_checkerboard(f(g/2) .+ f(3/2*g)') ./2)) for g in g_values]))
data = data[:, sortperm(data[1,:])]
for i in 1:20
    plot!(p, g_values, data[:,i])
end
display(p)
savefig("2024_06_13_mathieu_N=2_ed.pdf")