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


mmax, Nsites, Nbonds, Nsweep, e_cutoff, 
		SVD_error, gstart, delta_g, Ng,
        mbond, pairs, evod, angle, Estrength, 
		Nstates, output_filename, parity_symmetry_type,
		inversion_symmetry_type = get_input_data("input_quick.yml"; default_filename="test")

###{tdvp_filename}.h5 will be where the data is stored (written to on the fly while propagating every 5th sweep by default)
###ToDo: we might want to parse the name for that file from input

#Calculate kinetic matrix and x operator#
Ttmp = kinetic(mmax)
Xtmp = Xoperator(mmax)
Ytmp = Yoperator(mmax)
Uptmp = Upoperator(mmax)
Downtmp = Downoperator(mmax)
mInverttmp = MInversionOperator(mmax)
SmallEProjtmp = SmallEnergyProjector(mmax; m=1)

use_parity_symmetry = (parity_symmetry_type == "even" || parity_symmetry_type == "odd")
if use_parity_symmetry
    symmetry=parity_symmetry
else
    symmetry=trivial_symmetry
end

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

#Determine number of interaction pairs per starting site#
Nsecond = zeros(Int64,(Nsites-1))
for i=1:Nsites-1
	if pairs == "nearest"
		Nsecond[i]=i+1
	elseif pairs == "allpairs"
		Nsecond[i]=Nsites
	end
end

include("operators.jl")
include("observer.jl")

if evod == "dvr"
    fac1 = 1.0
    fac2 = 1.0
else 
    fac1 = -1.0
    fac2 = 1.0im
end

sites = siteinds("PlaRotor",Nsites;dim=Nspec, conserve_parity=use_parity_symmetry, conserve_L=false)

Random.seed!(1234)
if parity_symmetry_type == "even"
	psi = randomMPS(sites,[1 for i in 1:Nsites]; linkdims=50)
elseif parity_symmetry_type == "odd"
	psi = randomMPS(sites,[[1 for i in 1:Nsites-1]..., 2*mmax]; linkdims=50)
else
	psi = randomMPS(sites; linkdims=50)
end
# psi = MPS(sites, [1 for i in 1:Nsites])


sweeps = Sweeps(Nsweep)
maxdim!(sweeps,10,10,10,10,10,10,10,10,10,10,10,10,10,10,20,20,20,20,20,20,20,30,30,30,30,30, 30,30,30,30,30,40)
setcutoff!(sweeps, e_cutoff)

g = gstart
H = create_Hamiltonian(g, sites, Nsecond)
energy_eigenstates = MPS[]


filename = format(Format(output_filename), g, Nsites, parity_symmetry_type)
println(filename)
h5open(filename, "w") do file
	write(file,"N", Nsites)
	write(file,"mmax", mmax)
	write(file,"g", g)
	write(file, "bond_dim", get_maxdims(sweeps))
	write(file, "parity", parity_symmetry_type)
end


# finding excited state with DMRG
for i in 1:Nstates
    energy, ψ = dmrg(H,energy_eigenstates, psi, sweeps;outputlevel=1, weight=30)
    push!(energy_eigenstates, ψ)
	println("Excitation: ", i)

	h5open(filename, "r+") do file
		write(file, string("energy_eigenstates/", i), energy_eigenstates[i])
	end
end

# using Printf

# for i in 1:Nstates
#     write(f, string("energy_eigenstates/", i), energy_eigenstates[i])
# end

# save(@sprintf("../output_data/DMRG_runs/DMRG_g=%0.2f.jld", g),"energy_eigenstates", energy_eigenstates,
#     "N", Nsites,"mmax", mmax, "bond_dim", get_maxdims(sweeps))

