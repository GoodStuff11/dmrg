using ITensors
using ITensorNetworks
using Graphs
using Observers
using Printf
# using ITensors.HDF5
using TupleTools
using JLD
import Random

# using PyPlot
# push!(LOAD_PATH,pwd())
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
		inversion_symmetry_type = get_input_data("input_quick.yml"; default_filename="psi0_N6_g")

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

g = 1.1

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


function find_mapping(mat)
    mapping = []
    rows = size(mat)[2]
    for i=1:rows
        index = argmax(abs.(mat[i:end,i]))
        if abs(mat[index, i]) < 0.3
            push!(mapping, nothing)
        else
            push!(mapping, index)
        end
    end
    return mapping
end


sweeps = Sweeps(30)
maxdim!(sweeps,20)
setcutoff!(sweeps, 1e-10)

excitation_number = 20
Tij = zeros(excitation_number, excitation_number)
previous_eigenstates = MPS[]
prev_x = nothing
energies = []
g_values = []
connection = []


Hij = zeros(excitation_number, excitation_number)
Sij = zeros(excitation_number, excitation_number) # we want H x = lambda S x
tmp = zeros(excitation_number, excitation_number)

for g=0.1:0.05:1.6
    H = create_Hamiltonian(g, sites, Nsecond)
    energy_eigenstates = MPS[]

    push!(g_values, g)
    # finding excited state with DMRG
    for i in 1:excitation_number
        energy, ψ = dmrg(H,energy_eigenstates, psi, sweeps;outputlevel=0, weight=30)
        push!(energy_eigenstates, ψ)
    end
    # compute <i|H|j> and <i|j>
    for i in 1:excitation_number
        for j in 1:excitation_number
            Hij[i,j] = inner(energy_eigenstates[i], apply(H, energy_eigenstates[j]))
            Sij[i,j] = inner(energy_eigenstates[i], energy_eigenstates[j])
        end
    end
    # compute generalized eigenvalue problem to rediagonalize
    F = eigen(Hij, Sij)
    push!(energies, F.values)

    # comparing new eigenstates with previous to see which they map to
    if length(previous_eigenstates) > 0
        for i in 1:excitation_number
            for j in 1:excitation_number
                Tij[i,j] = inner(energy_eigenstates[i], previous_eigenstates[j])
            end
        end
        tmp = F.vectors

        push!(connection, find_mapping(F.vectors' * Tij *prev_x))
    else
        push!(connection, nothing)
    end
    prev_x = F.vectors
    previous_eigenstates = copy(energy_eigenstates)
end
connection = connection[2:end]
data = hcat(energies...)'


save("../output_data/DMRG_data.jld", "energies", data, "connections", 
    connection, "N", Nsites,"mmax", mmax, "bond_dim", get_maxdims(sweeps))