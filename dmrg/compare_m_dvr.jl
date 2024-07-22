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

push!(LOAD_PATH,pwd())
include("input_data.jl")
include("matrices.jl")
include("states.jl")
include("expectations.jl")
include("dvr.jl")
include("utility_funcs.jl")

Nspec, Nsites, Nbonds, Nsweep, e_cutoff, 
		SVD_error, gstart, delta_g, Ng,
        mbond, pairs, evod, angle, Estrength, 
		Nstates, output_filename, parity_symmetry_type,
		inversion_symmetry_type = get_input_data("input_quick_DMRG.yml"; default_filename="test")

Ttmp = kinetic(Nspec)
Xtmp = Xoperator(Nspec)
Ytmp = Yoperator(Nspec)
Uptmp = Upoperator(Nspec)
Downtmp = Downoperator(Nspec)

use_inversion_symmetry = inversion_symmetry_type == "even" || inversion_symmetry_type == "odd"
use_parity_symmetry = parity_symmetry_type == "even" || parity_symmetry_type == "odd"

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
end

include("operators.jl")
include("observer.jl")

sites = siteinds("PlaRotor",Nsites;basis=evod, dim=Nspec, conserve_parity=use_parity_symmetry, conserve_inversion_symmetry=use_inversion_symmetry)
    
function get_energy(gstart, evod, _Nspec;Nsweep=10, parity_symmetry_type = "even", inversion_symmetry_type = "even", Nsites = 3)
    ###{tdvp_filename}.h5 will be where the data is stored (written to on the fly while propagating every 5th sweep by default)
    ###ToDo: we might want to parse the name for that file from input

    #Calculate kinetic matrix and x operator#
    global Nspec = _Nspec
    
    Ttmp = kinetic(Nspec)
    Xtmp = Xoperator(Nspec)
    Ytmp = Yoperator(Nspec)
    Uptmp = Upoperator(Nspec)
    Downtmp = Downoperator(Nspec)

    use_inversion_symmetry = inversion_symmetry_type == "even" || inversion_symmetry_type == "odd"
    use_parity_symmetry = parity_symmetry_type == "even" || parity_symmetry_type == "odd"

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
    end

    include("operators.jl")
    include("observer.jl")
    
    sites = siteinds("PlaRotor",Nsites;basis=evod, dim=Nspec, conserve_parity=use_parity_symmetry, conserve_inversion_symmetry=use_inversion_symmetry)

    Random.seed!(1234)
    psi = generate_initial_state(sites; parity_symmetry_type, inversion_symmetry_type)
    # psi = MPS(sites, [1 for i in 1:Nsites])


    sweeps = Sweeps(Nsweep)
    maxdim!(sweeps,50)
    setcutoff!(sweeps, e_cutoff)

    g = gstart
    H = create_Hamiltonian(g, sites, "nearest"; evod=evod)

    energy, ψ = dmrg(H, psi, sweeps;outputlevel=0)
    dE = sqrt(abs(inner(ψ, apply(H, apply(H, ψ))) - energy^2))
    return energy, dE
end

using DataFrames, CSV
parity_symmetry_type = "none"
inversion_symmetry_type = "none"
nodenames = ["g", "basis", "Nspec", "parity", "inversion", "Nsweep", "Nsites", "energy", "energy uncertainty"]
df = DataFrame([[] for _ = nodenames] , nodenames)
(parity_symmetry_type,inversion_symmetry_type) = ("none","none")
for Nspec in reverse([7,8,9,10,11])
    for g in reverse([0.2,0.5,0.8,1.3,1.7,2])
        for Nsites in reverse([3,5,10])
            Nsweep =50
            for evod in ["m", "dvr"]
                try
                    @time energy = get_energy(g, evod, Nspec;Nsweep,parity_symmetry_type,inversion_symmetry_type, Nsites)
                    push!(df, [g, evod, Nspec, parity_symmetry_type, inversion_symmetry_type, Nsweep, Nsites, energy...])
                    @show g, evod, Nspec,Nsweep,parity_symmetry_type,inversion_symmetry_type, Nsites
                catch
                    continue
                end
            end
        end
    end
end

CSV.write("outputfile5.csv",df)
