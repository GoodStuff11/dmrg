using ITensors
# using ITensorNetworks
using Observers
using Printf: Format, format
using Printf
using TupleTools
using JLD
using YAML
using HDF5
import Random
using IterTools

push!(LOAD_PATH,pwd())
include("input_data.jl")
include("matrices.jl")
include("states.jl")
include("expectations.jl")
include("dvr.jl")
include("utility_funcs.jl")


path = "."
files = readdir(path)
even_files = sort(filter(x->occursin("even", x), files),by=x->split(x, "_"))
odd_files = filter(x->occursin("odd", x), files)

function apply_elementwise_projector(matrix, sites, state)
    tmp = copy(state)
    for i in eachindex(sites)
        tmp = apply(op(matrix,sites[i]),state)
    end
    return inner(state, (tmp + state)/2)
end
function apply_pairs_projector(matrix, sites, state)
    tmp = copy(state)
    for i in eachindex(sites)
        tmp = apply(op(matrix,sites[i], sites[length(sites)+1-i]),state)
    end
    return inner(state, (tmp + state)/2)
end


function compare_states(path, filename1, filename2; parity_symmetry="even")
    path1 = joinpath(path, filename1)
    path2 = joinpath(path, filename2)
    h5open(path1, "r") do f1
        # write(file, string("energy_eigenstates/", i), energy_eigenstates[i])
        Nsites = read(f1, "N")
        mmax = read(f1, "mmax")
    end

    use_parity_symmetry = (parity_symmetry_type == "even" || parity_symmetry_type == "odd")
    if use_parity_symmetry
        symmetry=parity_symmetry
    else
        symmetry=trivial_symmetry
    end

    #Define basis#
    mInvert = symmetry(MInversionOperator(mmax))
    refop = ReflectionOperator(mmax)

    use_parity_symmetry = (parity_symmetry_type == "even" || parity_symmetry_type == "odd")
    if use_parity_symmetry
        symmetry=parity_symmetry
    else
        symmetry=trivial_symmetry
    end

    # reading the mps to an array
    mps1_list = []
    mps2_list = []
    h5open(path1, "r") do f1        
        h5open(path2, "r") do f2
            for i = 1:20
                try
                    mps = read(f, @sprintf("energy_eigenstates/%d", i), MPS)
                    push!(mps1_list, mps)
                catch
                end
                try
                    mps = read(f, @sprintf("energy_eigenstates/%d", i), MPS)
                    push!(mps2_list, mps)
                end
            end
        end
    end

    overlap = zeros(length(mps1_list), length(mps2_list))
    for (i, mps1) in enumerate(mps1_list)
        for (j,mps2) in enumerate(mps2_list)
            overlap[i,j] = inner(mps1, mps2)
        end
    end
    println(overlap)
end

curves = [] # vector of vector with two elements. First is g values, second is y value
for (filename1, filename2) in partition(even_files,2,1)
    compare_states(path, filename1, filename2; parity_symmetry="even")
end