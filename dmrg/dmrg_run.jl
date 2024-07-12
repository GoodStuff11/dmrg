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
		inversion_symmetry_type = get_input_data("input_quick.yml"; default_filename="test")

###{tdvp_filename}.h5 will be where the data is stored (written to on the fly while propagating every 5th sweep by default)
###ToDo: we might want to parse the name for that file from input

#Calculate kinetic matrix and x operator#
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


Nspec=size(T,1)

include("operators.jl")
include("observer.jl")

sites = siteinds("PlaRotor",Nsites;basis=evod,dim=Nspec, conserve_parity=use_parity_symmetry, conserve_inversion_symmetry=use_inversion_symmetry)

Random.seed!(1234)
psi = generate_initial_state(sites; parity_symmetry_type, inversion_symmetry_type)


sweeps = Sweeps(Nsweep)
maxdim!(sweeps,10,10,10,10,10,10,10,10,10,10,10,10,10,10,20,20,20,20,20,20,20,30,30,30,30,30, 30,30,30,30,30,30,30, 30,35,35,35,35,35,35,35,35,35,50,50,50,50,50,50,50,50,50,50,60)
setcutoff!(sweeps, e_cutoff)

g = gstart
H = create_Hamiltonian(g, sites, "nearest"; evod=evod)
energy_eigenstates = MPS[]


filename = format(Format(output_filename), g, Nsites, parity_symmetry_type, inversion_symmetry_type)
println(filename)
h5open(filename, "w") do file
	write(file,"N", Nsites)
	write(file,"Nspec", Nspec)
	write(file,"g", g)
	write(file, "bond_dim", get_maxdims(sweeps))
	write(file, "parity", parity_symmetry_type)
	write(file, "inversion", inversion_symmetry_type)
	write(file, "basis", evod)
end


# finding excited state with DMRG
ground_energy = nothing
for i in 1:Nstates
    energy, ψ = dmrg(H,energy_eigenstates, psi, sweeps;outputlevel=1, weight=30)
    push!(energy_eigenstates, ψ)
	if i == 1
		global ground_energy = energy
	end
	println("Excitation: ", i)

	h5open(filename, "r+") do file
		write(file, string("energy_eigenstates/", i), energy_eigenstates[i])
		write(file, "energy", real.(energy))
	end
	if energy > ground_energy + 2*(Nsites-1)
		break
	end
end

# using Printf

# for i in 1:Nstates
#     write(f, string("energy_eigenstates/", i), energy_eigenstates[i])
# end

# save(@sprintf("../output_data/DMRG_runs/DMRG_g=%0.2f.jld", g),"energy_eigenstates", energy_eigenstates,
#     "N", Nsites,"mmax", mmax, "bond_dim", get_maxdims(sweeps))

