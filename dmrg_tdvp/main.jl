using ITensors
using ITensorTDVP
using Observers
using Printf
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

#Define output files#
create_file("energy.txt")
create_file("entropy_vN.txt")
create_file("entropy_Renyi.txt")
create_file("mux.txt")
create_file("muy.txt")
create_file("xcorr.txt")
create_file("ycorr.txt")
create_file("corr.txt")
for i=1:Nstates+1
	create_file("schmidt_values_"*string(i)*".txt")
end

for ig = 0:Ng-1
let
	include("operators.jl")
	include("observer.jl")

	if evod == "dvr"
		fac1 = 1.0
		fac2 = 1.0
	else 
		fac1 = -1.0
		fac2 = 1.0im
	end

	#Non-interacting rotors as initial guess#
	if ig == 0
		sites=siteinds("PlaRotor",Nsites)
		ampo0 = AutoMPO()
		#Define Hamiltonian as MPO#
		for j=1:Nsites
			ampo0 += 1.0,"T",j
			#Electric field#
			ampo0 += -cos(angle)*Estrength,"X",j
			ampo0 += -sin(angle)*Estrength*fac2,"Y",j
		end	
		H0 = MPO(ampo0,sites)
		#Define accuracy parameters#
		sweeps = Sweeps(Nsweep)
		#Set up initial state#
		global psi0 = randomMPS(sites,10)
		maxdim!(sweeps,10) # gradually increase states kept
		cutoff!(sweeps,SVD_error) # desired truncation error
		
		#Perform DMRG runs#
		obs = DemoObserver(e_cutoff)
		energy,psi = dmrg(H0,psi0,sweeps,observer=obs, outputlevel=0)
		
		global psi0 = psi
		global psit0 = copy(psi)

		maxbond=maxlinkdim(psi)
		f=open("log","a")
		println(f,"Max. bond dimension: ",maxbond)
		println(f)
		@printf(f,"Final energy = %.8f \n",energy)
		println(f)
		println(f,"Initial state calculated")
		println(f,"###############################################################################")
		println(f)
		close(f)
	end	
	g = gstart + ig*delta_g
	f=open("log","a")
	println(f,"##################################")
	println(f,"########g= ",g," ########")
	println(f,"##################################")
	println(f)
	println(f,"####DMRG calculation####")
	println(f,"Construct MPO")
	

	sites = siteinds(psi0)
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
			ampo += 1.0*c*fac1,"Y",i,"Y",j
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
	#Set up initial state#
	cutoff!(sweeps,SVD_error) # desired truncation error
	if ig == 0
		maxdim!(sweeps,10,20,30,40,Nbonds) # gradually increase states kept
	else	
		maxdim!(sweeps,maxlinkdim(psi0),Nbonds)
	end

    #Perform DMRG runs#
    println(f,"Start DMRG run")
    close(f)
    obs = DemoObserver(e_cutoff)
    energy,psi = dmrg(H,psi0,sweeps,observer=obs, outputlevel=0)

	psit0=apply(muX,psi)

	println("g= ",g," <psit0|H|psit0>=  ",inner(psit0,  Apply(H, psit0)) / inner(psit0, psit0))

	psi_ini=copy(psit0)
	time=0.
	ttotal=200
	dt=1.
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

	if g>0.1
	  
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

	psi0 = psi

	f=open("log","a")
	println(f)
	maxbond=maxlinkdim(psi)
	println(f,"Max. bond dimension: ",maxbond)
	println(f)
	@printf(f,"Final energy = %.8f \n",energy)
	println(f)
	close(f)

	if Nstates != 0

		energies = []
                append!(energies,energy)
                wavefunction = [psi for ii=1:Nstates+1]
                if ig == 0
                        global initial_states = [psi0 for ii=1:Nstates]
                end

		for istates=1:Nstates
			if ig == 0
	                	global initial_states[istates] = randomMPS(sites,Nbonds)
	        	else
	                	maxdim!(sweeps,maxlinkdim(initial_states[istates]),Nbonds) # gradually increase states kept
	        end
			cutoff!(sweeps,SVD_error) # desired truncation error
			energy ,psi = dmrg(H,wavefunction[1:istates],initial_states[istates] ,sweeps, observer=obs, outputlevel=0)
			global initial_states[istates] = psi
			f2=open("log","a")
			maxbond=maxlinkdim(psi)
			println(f,"Max. bond dimension: ",maxbond)
			println(f2)
			println(f2,"Final energy "*string(istates)*". excited state= "*string(round(energy,digits=12))*"\n")
			println(f2)
			close(f2)
			wavefunction[istates+1]=psi
			append!(energies,energy)
		end

		text_energy=" "
		text_ent_vN=" "
		text_ent_R=" "
		text_xcorr=" "
		text_ycorr=" "
		text_corr=" "
		text_mux=" "
		text_muy=" "
		for istates=1:Nstates+1
			
			#Calculate von-Neumann entropy and Schmidt coefficients#
			SvN,Renyi,Svalues = vN_entropy(wavefunction[istates],mbond)
		
			write_output("schmidt_values_"*string(istates)*".txt",g,Svalues)
			
			#Calculate dipole correlations#
			xcorr,ycorr = correlation(wavefunction[istates],Nsites,Nspec,evod,X,Y)
			
			#Calculate summed dipole moment and fluctuation#
			MuX,MuY = polarization(wavefunction[istates],Nsites,Nspec,evod,X,Y)
			
			
			text_energy*=string(energies[istates]," ")
			text_ent_vN*=string(SvN," ")
			text_ent_R*=string(Renyi," ")
			text_xcorr*=string(xcorr," ")
			text_ycorr*=string(ycorr," ")
			text_corr*=string((xcorr+ycorr)/(Nsites-1)," ")
			text_mux*=string(MuX," ")
			text_muy*=string(MuY," ")

		end	

		write_text("energy.txt",g,text_energy)
		write_text("entropy.txt",g,text_ent_vN)
		write_text("xcorr.txt",g,text_xcorr)
		write_text("ycorr.txt",g,text_ycorr)
		write_text("corr.txt",g,text_corr)
		write_text("mux.txt",g,text_mux)
		write_text("muy.txt",g,text_muy)

	else
		#Write energy to file#
		f=open("energy.txt","a")
		println(f,round(g,digits=4)," ",energy)
		close(f)
	
		#Write entropy and Schmidt values to file#
		SvN,Renyi,Svalues = vN_entropy(psi,mbond)
	
		write_output("entropy_vN.txt",g,SvN)
		write_output("entropy_Renyi.txt",g,Renyi)
		write_output("schmidt_values.txt",g,Svalues)
		
		#Write nearest-neighbour correlation to file#
		xcorr,ycorr = correlation(psi,Nsites,Nspec,evod,X,Y)
		
		write_output("xcorr.txt",g,xcorr)
		write_output("ycorr.txt",g,ycorr)
		write_output("corr.txt",g,(xcorr+ycorr)/(Nsites-1))
		
		
		#Write polarization to file#
		MuX,MuY = polarization(psi,Nsites,Nspec,evod,X,Y)
	
		write_output("mux.txt",g,MuX)
		write_output("muy.txt",g,MuY)
	end

end
end
f=open("log","a")
println(f)
println(f,"Calculation finished.")
close(f)
