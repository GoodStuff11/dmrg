using ITensors
using ITensorNetworks
using Observers
using Printf
using Printf: Format, format
using HDF5
using LinearAlgebra
import Random
	
function compute_timestep(psi, H; Nstd=3, freq_multiple=8)
	if H isa ITensorNetworks.TTN
		Hpsi = apply(H, psi; init=psi, nsweeps=3)
	else
		Hpsi = apply(H, psi)
	end
	H_ex = inner(psi, Hpsi)
	# @show H_ex
	H2_ex = inner(Hpsi, Hpsi)
	# @show sqrt(H2_ex - H_ex^2)
	peak_energy = freq_multiple*(abs(H_ex) + Nstd*sqrt(H2_ex - H_ex^2))
	return 2π/peak_energy
end

function apply_elementwise_projector(matrix, sites, state; parity="even")
    tmp = copy(state)
    for i in eachindex(sites)
        tmp = apply(op(matrix,sites[i]),tmp)
    end
    if parity == "even"
        return (tmp + state)/2
    end
    return (state - tmp)/2
end


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
		inversion_symmetry_type = get_input_data("input_quick.yml"; default_filename="psi0_N6_g")

###{tdvp_filename}.h5 will be where the data is stored (written to on the fly while propagating every 5th sweep by default)
###ToDo: we might want to parse the name for that file from input

use_inversion_symmetry = inversion_symmetry_type == "even" || inversion_symmetry_type == "odd"
use_parity_symmetry = parity_symmetry_type == "even" || parity_symmetry_type == "odd"

#Calculate kinetic matrix and x operator#
Ttmp = kinetic(Nspec)
Xtmp = Xoperator(Nspec)
Ytmp = Yoperator(Nspec)
Uptmp = Upoperator(Nspec)
Downtmp = Downoperator(Nspec)

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


let
	starting_memory = total_memory_usage()
	include("operators.jl")
	include("observer.jl")

	g = 1.1

	tdvp_filename = format(Format(output_filename), g, Nsites, parity_symmetry_type, inversion_symmetry_type)
	# dmrg_filename = @sprintf("dmrg_logs/%sg=%.3f",output_filename,g)

	sites = siteinds("PlaRotor",Nsites;basis=evod,dim=Nspec, conserve_parity=use_parity_symmetry, conserve_inversion_symmetry=use_inversion_symmetry)

	Random.seed!(1234)
	psit0 = generate_initial_state(sites; parity_symmetry_type, inversion_symmetry_type)
	@show flux(psit0)
	# m_inv = MPO(create_transform_ampo("mInvert"), sites)
	# even_projection = MPO(create_transform_ampo("mInvert"; coefficient=0.5)+(0.5,"Id",1), sites)
	# odd_projection = MPO(create_transform_ampo("mInvert"; coefficient=-0.5)+(0.5,"Id",1),sites)
	# psit0 = apply_elementwise_projector(symmetry(MParityOperator(Nspec)), sites, psi; parity=inversion_symmetry_type)
	# normalize!(psit0)
	# @show typeof(psit0)

	# sites = siteinds(psi)

	@show hasqns.(sites)
	H, ampo = create_Hamiltonian(g, sites, "nearest"; evod=evod, get_ampo=true)

	# println("g= ",g," before transform <psi|H|psi>=  ",inner(psi,  apply(H, psi)))
	println("g= ",g," after transform <psit0|H|psit0>=  ",inner(psit0',  apply(H, psit0)))

	#################################
	# prepare TTN of state
	psi_ini = MPS_to_ITensorNetwork(psit0)
	Httn = MPS_to_ITensorNetwork(H)
	SvN, Renyi, Svalues = vN_entropy(psit0)
	println(SvN, " ", Svalues)
	# fct=open(string("tdvp_logs/itn_ct",string(round(g,digits=3))),"w")
	
	
	# prepare exponential decay state
	nsweeps = 4
	dt = round(20*compute_timestep(psit0, H; freq_multiple=4), sigdigits=5)
	ttotal = round(dt*nsweeps,digits=8)
	dt = ttotal/nsweeps  + 1e-15 # 1e-5 prevents the ceil(ttotal/dt) from being not equal to nsweeps
	@show ceil(real(ttotal/dt))
	@show dt
	
	
	# @show expect(ampo, ITensorNetworks(psi_ini))
	@show inner(psit0', H, psit0)
	
	psi_ini = ITensorNetworks.tdvp(
		ITensorNetworks.TTN(Httn),
		-ttotal,
		ITensorNetworks.TTN(psi_ini);
		# time_step=-dt,
		nsweeps=nsweeps,
		outputlevel=1,
		normalize=true,
		maxdim=50,
		cutoff=1e-6,
		nsites=1, # 2 alows maxlinkdim to increase
	)
	@show inner(TTN_to_MPS(psi_ini)', H, TTN_to_MPS(psi_ini))
	@show vN_entropy(TTN_to_MPS(psi_ini))[3]

	##################################
	# applying tdvp

	# adjusting the Hamiltonian by a constant factor to have a range of energies 
	# around zero
	H_ex = inner(TTN_to_MPS(psi_ini)', H, TTN_to_MPS(psi_ini))
	tmp = MPO(ampo + (-H_ex, "Id",1),sites)
	new_Httn = ITensorNetworks.TTN(MPS_to_ITensorNetwork(tmp))
	
	nsweeps = Nsweep
	dt = round(compute_timestep(psi_ini, new_Httn; freq_multiple=8, Nstd=2), sigdigits=10)
	ttotal = round(dt*nsweeps,digits=8)
	dt = ttotal/nsweeps  + 1e-15
	@show dt

	
	function get_memory()
		return total_memory_usage() - starting_memory
	end

	function step(; which_sweep)
		return which_sweep
	end
	
	function return_time(; which_sweep)
		return abs.(which_sweep*dt)
	end
	
	  
	function step_duration(; sweep_time)
		return sweep_time
	end
	
	function return_Svals(;state)
		return vN_entropy(TTN_to_MPS(state))[3]
	end
	function return_maxlinkdim(;state)
		return maxlinkdim(state)
	end 
	
	function return_corr(;state)
		return correlation(TTN_to_MPS(state), evod)
	end
	function return_overlap(; state)
		return inner(psi_ini, state)
	end
	
	obs = observer(
		"steps" => step, "time" => return_time,"step_duration"=>step_duration, 
		"overlap" => return_overlap, "corr" => return_corr, "Svals" => return_Svals, "maxlinkdim" => return_maxlinkdim,
		"memory(GB)" => get_memory
		# "norm" => return_norm
	)

	@show maxlinkdim(psi_ini) 


	psi_f = ITensorNetworks.tdvp(
		new_Httn,
		-im * ttotal,
		psi_ini;
		time_step=-im * dt,
		outputlevel=1,
		normalize=true,
		maxdim=Nbonds,
		cutoff=e_cutoff,
		nsites=2,
		(sweep_observer!)=obs,
	)
	println(obs)
	savedata(tdvp_filename, obs; H_ex=H_ex)

end
