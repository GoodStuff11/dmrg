using ITensors
using ITensorTDVP
using Observers
using Printf
using HDF5
push!(LOAD_PATH,pwd())
using matrices
using states
using expectations
using dvr
using utility_funcs

import YAML

data = YAML.load_file("input_quick.yml")
mmax = data["m"]
Nsites = data["Nsites"]
Nbonds = data["Nbond"]
Nsweep = data["Nsweep"]
e_cutoff = data["ecutoff"]
SVD_error = data["SVD"]
gstart = data["gstart"]
delta_g = data["dg"]
Ng = data["Ng"]
mbond = data["mbond"]
pairs = data["pairs"]
evod = data["states"]
angle = data["angle"]*pi/180.0
Estrength = data["strength"]
Nstates = data["Nstates"]
output_path = data["outputpath"]

f=open("log","w")

println(f,"#######################################")
println(f,"###########Basis information###########")
println(f,"#######################################")
println(f,"mmax= ",mmax)
println(f,"Number of sites: ",Nsites)
println(f,"Dimension of local Hilbert space: ",2*mmax+1)

#Calculate kinetic matrix and x operator#
Ttmp = kinetic(mmax)
Xtmp = Xoperator(mmax)
Ytmp = Yoperator(mmax)

#Define basis#
if evod == "all"
	global T = Ttmp
	global X = Xtmp
	global Y = Ytmp
	Nspec=size(T,1)
	println(f,"all m-states are considered")
elseif evod == "dvr"
	tmp1,tmp2,tmp3 = exp_dvr(mmax)
	global T = tmp1
	global X = tmp2
	global Y = tmp3
	Nspec=size(T,1)
	println(f,"DVR-basis is used")
else
	global T = symmetry(Ttmp,evod,mmax)
	global X = symmetry(Xtmp,evod,mmax)
	global Y = symmetry(Ytmp,evod,mmax)
	println(f,"only "+evod+" m-states are considered")
	Nspec=size(T,1)
end
println(f,"Dimension of local Hilbert space for chosen m-states: ",Nspec)
if pairs == "nearest"
	println(f,"only nearest-neighbour interactions")
elseif pairs == "allpairs"
	println(f,"all interactions")
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

println(f)
println(f,"#################################################################################")
println(f,"##################################")
println(f,"####Calculate free rotor chain####")
println(f,"##################################")
println(f)
close(f)

let

	include("dmrg_tdvp/operators.jl")
	include("dmrg_tdvp/observer.jl")

	if evod == "dvr"
		fac1 = 1.0
		fac2 = 1.0
	else 
		fac1 = -1.0
		fac2 = 1.0im
	end

	g = 1.1

	f=open("log","a")
	println(f,"##################################")
	println(f,"########g= ",g," ########")
	println(f,"##################################")
	println(f)
	println(f,"####DMRG calculation####")
	println(f,"Construct MPO")
	
	sites=siteinds("PlaRotor",Nsites)

	mps_out=h5open(string("psi0_g",string(round(g,digits=3))),"r")
	psi=read(mps_out,"MPS",MPS)
	close(mps_out)

	sites = siteinds(psi)

	ampo = AutoMPO() 
	dipoleX = AutoMPO() 
	dipoleY = AutoMPO() 
	for i=1:Nsites-1
		ampo += 1.0,"T",i
		dipoleX+=1.0,"X",i
		dipoleY+=1.0,"Y",i
		for j=i+1:Nsecond[i]
			c=g/((abs(j-i))^3)
			#y_iy_j#
			ampo += -1.0*c*fac1,"Y",i,"Y",j
			#2*x_ix_j#
			ampo += -2.0*c,"X",i,"X",j
		end
		#Electric field#
		ampo += -cos(angle)*Estrength,"X",i
		ampo += -sin(angle)*Estrength*fac2,"Y",i
	end
	ampo += 1.0,"T",Nsites
	dipoleX+=1.0,"X",Nsites
	dipoleY+=1.0,"Y",Nsites
	#Electric field#
	ampo += -cos(angle)*Estrength,"X",Nsites
	ampo += -sin(angle)*Estrength*fac2,"Y",Nsites

	H = MPO(ampo,sites)
	muX = MPO(dipoleX,sites)
	muY = MPO(dipoleY,sites)
	#Define accuracy parameters#
	sweeps = Sweeps(Nsweep)

    #Perform TDVP

	psit0=apply(muX,psi)

	println("g= ",g," <psit0|H|psit0>=  ",inner(psit0,  Apply(H, psit0)) / inner(psit0, psit0))

	psi_ini=copy(psit0)
	time=0.
	ttotal=200
	dt=.5
	fct=open(string("ct",string(round(g,digits=3))),"w")
	println(fct,time," ",1.," ",0.)

	function step(; sweep, bond, half_sweep)
		if bond == 1 && half_sweep == 2
		  return sweep
		end
		return nothing
	end
	  
	function current_time(; current_time, bond, half_sweep)
		if bond == 1 && half_sweep == 2
		  return current_time
		end
		return nothing
	end
	  
	function return_state(; psi, bond, half_sweep)
		if bond == 1 && half_sweep == 2
		  return psi
		end
		return nothing
	end
	function return_corr(; psi0, psi, bond, half_sweep)
		if bond == 1 && half_sweep == 2
		  return inner(psi, psi) / sqrt(inner(psi, psi)*inner(psi0, psi0))
		end
		return nothing
	end
	  
	obs2 = observer(
		"steps" => step, "times" => current_time, "psis" => return_state
	)

	
	psi_f = tdvp(
		H,
		-im * ttotal,
		psi_ini;
		time_step=-im * dt,
		outputlevel=0,
		normalize=false,
		maxdim=20,
		(observer!)=obs2,
	)
	
	steps = obs2.steps
	times = obs2.times
	psis = obs2.psis
	

	for n in 1:length(steps)	
		#psit = tdvp(H,-im*dt,psi_ini;nsweeps=10,reverse_step=false,normalize=true)
		#psit = tdvp(H,-im*dt,psi_ini;nsweeps=10,reverse_step=false,normalize=true,maxdim=50,cutoff=1e-10)
		#psit = tdvp(H,-1.0,psi_ini;nsweeps=10,reverse_step=false,normalize=true)
		#time=time+dt
		#tcorr=inner(psit, psit0)
		tcorr=inner(psis[n], psit0) / sqrt(inner(psis[n], psis[n])*inner(psit0, psit0))
		println(fct,-imag(times[n])," ",real(tcorr)," ",imag(tcorr))
		#psi_ini=copy(psit)
	end
	println(fct,"\n")
	close(fct)


end

f=open("log","a")
println(f)
println(f,"Calculation finished.")
close(f)
