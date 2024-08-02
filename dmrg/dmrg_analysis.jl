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
        tmp = apply(op(matrix,sites[i]),tmp)
    end
    if parity == "even"
        return inner(state, (tmp + state)/2)
    end
    return inner(state, (tmp - state)/2)
end
function apply_pairs_projector(matrix, sites, state; parity="even")
    tmp = copy(state)
    for i in 1:length(sites)÷2
        tmp = apply(op(matrix,sites[i], sites[end+1-i]),tmp)
    end
    if parity == "even"
        return inner(state, (tmp + state)/2)
    end
    return inner(state, (state - tmp)/2)
end

function get_mps_list(path)
    mps_list = []
    energy = []
    iteration_enegy = []
    energy_uncertainty = []
    h5open(path, "r") do file   
        for i = 1:20
            println(i)
            evod = read(file, "basis")
            Nspec = read(file, "Nspec")
            try
                mps = read(file, "energy_eigenstates/$i", MPS)
                push!(energy,read(file, "energy/$i"))
                push!(iteration_enegy,read(file, "Delta H/$i"))
                push!(energy_uncertainty,read(file, "iteration_energy/$i"))
                push!(mps_list, mps)
            catch LoadError
                println("That's all folks")
                break
            end
        end
    end    
    return mps_list, energy, iteration_enegy, energy_uncertainty
end


function compute_states(path, filename, past_mps_list, past_vectors; parity_symmetry_type="even",inversion_symmetry_type="even", dim=11, evod="m")
    path = joinpath(path, filename)

    Ttmp = kinetic(dim)
    Xtmp = Xoperator(dim)
    Ytmp = Yoperator(dim)
    Uptmp = Upoperator(dim)
    Downtmp = Downoperator(dim)
    
    use_inversion_symmetry = inversion_symmetry_type == "even" || inversion_symmetry_type == "odd"
    use_parity_symmetry = parity_symmetry_type == "even" || parity_symmetry_type == "odd"

    #Define basis
    if evod == "dvr"
        symmetry = trivial_symmetry
        if use_inversion_symmetry && use_parity_symmetry
            symmetry = dvr_symmetric_basis
        elseif use_inversion_symmetry
            symmetry = dvr_inversion_symmetry
        elseif use_parity_symmetry
            symmetry = dvr_rotation_symmetry
        end
    
        # define basis
        tmp1,tmp2,tmp3 = symmetry.(exp_dvr(Nspec))
        global T = tmp1
        global X = tmp2
        global Y = tmp3

        mInvert = symmetry(phiReflectionOperator(dim))
    end
    if evod == "m"
        symmetry = trivial_symmetry
        if use_inversion_symmetry && use_parity_symmetry
            symmetry = x -> parity_symmetry(m_inversion_symmetry(x))
        elseif use_inversion_symmetry
            symmetry = m_inversion_symmetry
        elseif use_parity_symmetry
            symmetry = parity_symmetry
        end
    
        # define basis
        global T = symmetry(Ttmp)
        global X = symmetry(Xtmp)
        global Y = symmetry(Ytmp)
        global Up = symmetry(Uptmp)
        global Down = symmetry(Downtmp)

        mInvert = symmetry(MInversionOperator(dim))
    end
    

    #Define basis#
    
    # refop = ReflectionOperator(dim)

    # reading the mps to an array
    println("Getting MPS")
    mps_list, energy, iteration_enegy, energy_uncertainty = get_mps_list(path)

    # for (i, mps) in enumerate(mps_list)
    #     println(i)
    #     println(apply_elementwise_projector(mInvert, sites, mps))
    #     println(apply_pairs_projector(refop, sites, mps))
    # end

    g = parse(Float64, split(filename, "_")[2][3:end])
    sites = siteinds(mps_list[1])
    H = create_Hamiltonian(g, sites, "nearest")

    Hij = zeros((length(mps_list), length(mps_list)))
    Sij = zeros((length(mps_list), length(mps_list))) # we want H x = lambda S x

    println("Computing H and S")
    for (i, mps1) in enumerate(mps_list)
        for (j,mps2) in enumerate(mps_list)
            Sij[i,j] = inner(mps1, mps2)
            Hij[i,j] = inner(mps1, apply(H,mps2))
        end
    end
    println("Computing eigenvalues")
    F = eigen(Hij, Sij)
    energy_levels = F.values
    # eigen_vectors = F.vectors ./ sqrt.(diag(F.vectors' * F.vectors)') # normalize

    # println("Computing overlap")
    # overlap = nothing
    # if !isnothing(past_mps_list)
    #     overlap = zeros((length(mps_list), length(past_mps_list)))
    #     for (i, mps1) in enumerate(mps_list)
    #         for (j,mps2) in enumerate(past_mps_list)
    #             overlap[i,j] = inner(mps1, mps2)
    #         end
    #     end

    #     overlap = eigen_vectors' * overlap * past_vectors # adjustment may not be necessary depending on DMRG accuracy
    # end

    # println("Computing parity")
    # even_m_parity = [apply_elementwise_projector(mInvert, sites, mps; parity="even") for mps in mps_list]
    # even_reflection_parity = [apply_pairs_projector(refop, sites, mps) for mps in mps_list] # too computationally extensive
    iteration_uncertainty = iteration_enegy[end] - iteration_enegy[end - 1]
    return g, energy_levels, iteration_uncertainty, energy_uncertainty, Sij, Hij
end

function group_files_in_dir(path)
    # groups all files in folder into their symmetries and sorts them by g value
    # requires the naming convention 
    # SOMETHING_g=#_N=#_parity=even_inversion=odd.jld 
    files = readdir(path)
    filename_groupings = Dict()
    for file in files
        parity_key = nothing
        inversion_key = nothing
        if occursin("parity=even", file)
            parity_key = "even"
        elseif occursin("parity=odd", file)
            parity_key = "odd"
        else
            parity_key = "none"
        end
    
        if occursin("inversion=even", file)
            inversion_key = "even"
        elseif occursin("inversion=odd", file)
            inversion_key = "odd"
        else
            inversion_key = "none"
        end
    
        if (parity_key, inversion_key) in keys(filename_groupings)
            push!(filename_groupings[(parity_key, inversion_key)], file)
        else
            filename_groupings[(parity_key, inversion_key)] = [file]
        end
    end
    for (key, val) in filename_groupings
        filename_groupings[key] = sort(val, by=x->parse(Float64, split(x, "_")[2][3:end]))
    end
    return filename_groupings
end

include("operators.jl")
include("observer.jl")

_, _, _, _, _, 
		_, gstart, _, _,
        _, _, _, _, _, 
		_, _, parity_symmetry_type,
		inversion_symmetry_type = get_input_data("input_quick_DMRG.yml"; default_filename="test")

path = raw"/home/jkambulo/projects/def-pnroy/jkambulo/dmrg/output_data/DMRG_runs3/"
# path = raw"C:\Users\jonat\OneDrive\Documents\programming\AnacondaProjects\PHYS437B\dmrg\output_data\DMRG_runs_test"
filenames = readdir(path)

Nspec = nothing
Nsites = nothing
evod = nothing
overlap_threshold = 0.4
curves = Dict() # vector of vector with two elements. First is g values, second is y value
_past_mps_list = nothing
_past_vectors = nothing
num_completed_curves = 0
println(path)
for filename in filenames

    h5open(joinpath(path, files[1]), "r") do f1
        # write(file, string("energy_eigenstates/", i), energy_eigenstates[i])
        global Nsites = read(f1, "N")
        global Nspec = read(f1, "Nspec")
        global evod = read(f1, "basis")
        global parity_symmetry_type = read(file, "parity")
        global inversion_symmetry_type = read(file, "inversion")
    end

    ### delete me ###
    if !(occursin(@sprintf("g=%.2f", gstart), filename))
        continue
    end
    ####
    # dict_key(x) = "$parity_symmetry-$inversion_symmetry-$x"
    println(filename)
    g, energy_levels, iteration_uncertainty, 
        energy_uncertainty,Sij, Hij = compute_states(path, filename, _past_mps_list, _past_vectors; 
                                                parity_symmetry_type=parity_symmetry_type, 
                                                inversion_symmetry_type=inversion_symmetry_type, dim=Nspec, evod=evod)

    h5open(@sprintf("/home/jkambulo/projects/def-pnroy/jkambulo/dmrg/output_data/processed_data2/processed_data_%s", filename), "w") do file
        write(file, "g", g)  
        write(file, "energy_levels", energy_levels)  
        write(file, "iteration_uncertainty", iteration_uncertainty)
        write(file, "energy_uncertainty", energy_uncertainty)
        write(file, "Sij", Sij)
        write(file, "Hij", Hij)
    end
    

    # # adjusting curves based on overlap
    # if isnothing(overlap)
    #     for (i,(energy, parity)) in enumerate(zip(energy_levels, even_m_parity))
    #         curves[dict_key(i)] = [[g], [energy], [parity]]
    #     end
    # else
    #     overlap = abs.(overlap).^2
    #     # show(IOContext(stdout, :limit=>false), MIME"text/plain"(), overlap)
    #     new_size, old_size = size(overlap)
    #     used_indices = Set()
    #     mapping = Dict()
    #     for i=1:old_size
    #         max_index = argmax(overlap[:,i])
    #         if overlap[max_index, i] >= overlap_threshold && !(max_index in used_indices)
    #             mapping[dict_key(i)] = dict_key(max_index)
    #             push!(used_indices, dict_key(max_index))
    #         else
    #             mapping[dict_key(i)] = dict_key("curve$num_completed_curves")
    #             global num_completed_curves += 1
    #         end
    #     end
    #     global curves = Dict(get!(mapping, key, key)=>val for (key, val) in curves)
    #     for (i, (energy, parity)) in enumerate(zip(energy_levels, even_m_parity))
    #         if !(dict_key(i) in used_indices)
    #             curves[dict_key(i)] = [[g], [energy], [parity]]
    #         else
    #             push!(curves[dict_key(i)][1],g)
    #             push!(curves[dict_key(i)][2],energy)
    #             push!(curves[dict_key(i)][3],parity)
    #         end
    #     end

    # end
end
# println(curves)
# save(raw"/home/jkambulo/projects/def-pnroy/jkambulo/dmrg/output_data/processed_data/plot_data.jld", "curves", curves) 
