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
using LinearAlgebra

push!(LOAD_PATH,pwd())
include("input_data.jl")
include("matrices.jl")
include("states.jl")
include("expectations.jl")
include("dvr.jl")
include("utility_funcs.jl")

function apply_elementwise_projector(matrix, sites, state; parity="even")
    tmp = copy(state)
    for i in eachindex(sites)
        tmp = apply(op(matrix,sites[i]),state)
    end
    if parity == "even"
        return inner(state, (tmp + state)/2)
    end
    return inner(state, (tmp - state)/2)
end
function apply_pairs_projector(matrix, sites, state; parity="even")
    tmp = copy(state)
    for i in eachindex(sites)
        tmp = apply(op(matrix,sites[i], sites[length(sites)+1-i]),state)
    end
    if parity == "even"
        return inner(state, (tmp + state)/2)
    end
    return inner(state, (tmp - state)/2)
end

function get_mps_list(path)
    mps_list = []
    h5open(path, "r") do f1   
        for i = 1:20
            println(i)
            try
                mps = read(f1, @sprintf("energy_eigenstates/%d", i), MPS)
                push!(mps1_list, mps)
            catch LoadError
                break
            end
        end
    end    
    return mps_list 
end

function compute_states(path, filename, past_mps_list, past_vectors; parity_symmetry_type="even", mmax=5)
    path = joinpath(path, filename)
    
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
    mps_list = get_mps_list(path)

    # for (i, mps) in enumerate(mps_list)
    #     println(i)
    #     println(apply_elementwise_projector(mInvert, sites, mps))
    #     println(apply_pairs_projector(refop, sites, mps))
    # end

    g = parse(Float64, split(filename1, "_")[2][3:end])
    sites = siteinds(mps_list[1])
    H = create_Hamiltonian(g, sites, pairs)

    Hij = zeros((length(mps_list), length(past_mps_list)))
    Sij = zeros((length(mps_list), length(past_mps_list))) # we want H x = lambda S x
    tmp = zeros((length(mps_list), length(past_mps_list)))


    for (i, mps1) in enumerate(mps_list)
        for (j,mps2) in enumerate(mps_list)
            Sij[i,j] = inner(mps1, mps2)
            Hij[i,j] = inner(mps1, apply(H,mps2))
        end
    end
    F = eigen(Hij, Sij)
    energy_levels = F.values

    overlap = nothing
    if isnothing(past_mps_list)
        overlap = zeros((length(mps_list), length(past_mps_list)))
        for (i, mps1) in enumerate(mps_list)
            for (j,mps2) in enumerate(past_mps_list)
                overlap[i,j] = inner(mps1, mps2)
            end
        end
        overlap = F.vectors' * overlap * past_vectors # adjustment may not be necessary depending on DMRG accuracy
    end


    even_m_parity = [apply_elementwise_projector(mInvert, sites, mps; parity="even") for mps in mps_list]
    # even_reflection_parity = [apply_pairs_projector(refop, sites, mps) for mps in mps_list] # too computationally extensive

    return g, mps_list, adjusted_overlap, energy_levels, even_m_parity, F.vectors
end


path = raw"D:\datasets\DMRG_runs2"
files = readdir(path)
even_files = sort(filter(x->occursin("even", x), files),by=x->parse(Float64, split(x, "_")[2][3:end]))
odd_files = sort(filter(x->occursin("odd", x), files), by=x->parse(Float64, split(x, "_")[2][3:end]))


mmax = nothing
Nsites = nothing
h5open(joinpath(path, even_files[1]), "r") do f1
    # write(file, string("energy_eigenstates/", i), energy_eigenstates[i])
    global Nsites = read(f1, "N")
    global mmax = read(f1, "mmax")
end

include("operators.jl")
include("observer.jl")
Ttmp = kinetic(mmax)
Xtmp = Xoperator(mmax)
Ytmp = Yoperator(mmax)
Uptmp = Upoperator(mmax)
Downtmp = Downoperator(mmax)
mInverttmp = MInversionOperator(mmax)
SmallEProjtmp = SmallEnergyProjector(mmax; m=1)

evod = "m"
#Define basis#
if evod == "dvr"
	tmp1,tmp2,tmp3 = symmetry.(exp_dvr(mmax))
	global T = symmetry(tmp1)
	global X = symmetry(tmp2)
	global Y = symmetry(tmp3)

	Nspec=size(T,1)
else 
	global T = symmetry(Ttmp)
	global X = symmetry(Xtmp)
	global Y = symmetry(Ytmp)
	global Up = symmetry(Uptmp)
	global Down = symmetry(Downtmp)
	global mInvert = symmetry(mInverttmp)
	global SmallEProj = symmetry(SmallEProjtmp)

	Nspec=size(T,1)
end



# println(even_files)
# println(odd_files)
overlap_threshold = 0.4
curves = Dict() # vector of vector with two elements. First is g values, second is y value
past_mps_list = nothing
past_vectors = nothing
num_completed_curves = 0
for filename in even_files
    g, mps_list, overlap, energy_levels, even_m_parity, past_vectors = compute_states(path, filename, past_mps_list, past_vectors; parity_symmetry_type="even", mmax=mmax)

    # adjusting curves based on overlap
    if isnothing(overlap)
        for (energy, parity) in zip(energy_levels, even_m_parity)
            curves[i] = [[g], [energy], [parity]]
        end
    else
        new_size, old_size = size(overlap)
        used_indices = Set()
        mapping = Dict()
        for i=1:old_size
            max_index = argmax(overlap[:,i])
            if overlap[max_index, i] >= overlap_threshold && !(max_index in used_indices)
                mapping[i] = max_index
                push!(used_indices, max_index)
            else
                mapping[i] = string("curve", num_completed_curves)
                num_completed_curves += 1
            end
        end
        curves = Dict(mapping[key]=>val for (key, val) in curves)
        for (i, (energy, parity)) in enumerate(zip(energy_levels, even_m_parity))
            if !(i in used_indices)
                curves[i] = [[g], [energy], [parity]]
            else
                push!(curves[i][1],g)
                push!(curves[i][2],energy)
                push!(curves[i][3],parity)
            end
        end

    end
    break
end