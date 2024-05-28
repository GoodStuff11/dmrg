using ITensors
using Observers
using Printf
using ITensors.HDF5
using LinearAlgebra
#using StatsPlots
#push!(LOAD_PATH,pwd())
include("input_data.jl")
include("matrices.jl")
include("states.jl")
include("expectations.jl")
include("dvr.jl")
include("utility_funcs.jl")

mmax, Nsites, Nbonds, Nsweep, e_cutoff, 
		SVD_error, gstart, delta_g, Ng,
        mbond, pairs, evod, angle, Estrength, 
		Nstates, dmrg_output_filename = get_input_data("input_quick.yml"; 
											default_filename="psi0_N6_g")
@show mmax, Nsites, Nbonds, Nsweep, e_cutoff, 
	SVD_error, gstart, delta_g, Ng,
	mbond, pairs, evod, angle, Estrength, 
	Nstates, dmrg_output_filename
#f=open("log","w")
global f=stdout
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
Uptmp = Upoperator(mmax)
Downtmp = Downoperator(mmax)
#Define basis#
#matshow(Ttmp)
#matshow(Xtmp)
#show()
use_symmetry=false
if use_symmetry
    symmetry=parity_symmetry
else
    symmetry=trivial_symmetry
end
if evod == "all"
	global T = symmetry(Ttmp)
	global X = symmetry(Xtmp)
	global Y = symmetry(Ytmp)
	global Up = symmetry(Uptmp)
	global Down = symmetry(Downtmp)
	#matshow(T)
	#matshow(X)
	#show()
	#hilbert space size
	Nspec=size(T,1)
	println("all m-states are considered")
	println(f,"all m-states are considered")
elseif evod == "dvr"
	tmp1,tmp2,tmp3 = exp_dvr(mmax)
	global T = tmp1
	global X = tmp2
	global Y = tmp3
	Nspec=size(T,1)
	println("DVR-basis is used")
	println(f,"DVR-basis is used")
else
	global T = symmetry(Ttmp,evod,mmax)
	global X = symmetry(Xtmp,evod,mmax)
	global Y = symmetry(Ytmp,evod,mmax)
	println("only "+evod+" m-states are considered")
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
#close(f)

#Define output files#
create_file("dmrg_logs/energy.txt",evod,mmax,Nsites)
create_file("dmrg_logs/entropy_vN.txt",evod,mmax,Nsites)
create_file("dmrg_logs/entropy_Renyi.txt",evod,mmax,Nsites)
create_file("dmrg_logs/mux.txt",evod,mmax,Nsites)
create_file("dmrg_logs/muy.txt",evod,mmax,Nsites)
create_file("dmrg_logs/xcorr.txt",evod,mmax,Nsites)
create_file("dmrg_logs/ycorr.txt",evod,mmax,Nsites)
create_file("dmrg_logs/corr.txt",evod,mmax,Nsites)
f_int=open("dmrg_logs/int.txt","w")
for i=1:Nstates+1
	create_file("dmrg_logs/schmidt_values_"*string(i)*".txt",evod,mmax,Nsites)
end

for ig = 0:Ng-1
let
	f=stdout
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
		sites=siteinds("PlaRotor",Nsites;dim=Nspec, conserve_parity=use_symmetry, conserve_L=false)
		ampo0 = AutoMPO()
		#Define Hamiltonian as MPO#
		
		for j=1:Nsites
			ampo0 += 1.0,"T",j
			#Electric field#
			if !iszero(Estrength)
				ampo0 += -cos(angle)*Estrength,"X",j
				ampo0 += -sin(angle)*:qngth*fac2,"Y",j
			end
		end	
		@show "created opsum"
		H0 = MPO(ampo0,sites)
		
		@show "created MPO"
		#Define accuracy parameters#
		sweeps = Sweeps(Nsweep)
		#Set up initial state#
		target_sector=0
		v=[zeros(Nspec) for i in 1:length(sites)]
		for i in eachindex(v)
			v[i][rand(1:Nspec)]=1.0
		end
		# psi = randomMPS(sites,(i)->rand(1:2:Nspec); linkdims=10)
		# mps_out=h5open("random_N6","w")
		# write(mps_out,"MPS",psi)
		# close(mps_out)
		# throw("ERROE")
	
		global psi0 = randomMPS(sites,(i)->rand(1:2:Nspec))
		
		maxdim!(sweeps,10) # gradually increase states kept
		cutoff!(sweeps,SVD_error) # desired truncation error
		
		#Perform DMRG runs#
		obs = DemoObserver(e_cutoff)
		energy,psi = dmrg(H0,psi0,sweeps,observer=obs, outputlevel=0)
		@show energy
		global psi0 = psi

		maxbond=maxlinkdim(psi)
		f=open("dmrg_logs/log","a")
		println(f,"Max. bond dimension: ",maxbond)
		println(f)
		@printf(f,"Final energy = %.8f \n",energy)
		println(f)
		println(f,"Initial state calculated")
		println(f,"###############################################################################")
		println(f)
		#close(f)
	end
	@show("done with trivial calc")
	g = gstart + ig*delta_g
	#f=open("log","a")
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
			#ampo += 1.0*c*fac1,"Y",i,"Y",j
			#2*x_ix_j#
			#ampo += -2.0*c,"X",i,"X",j
			# up up
			ampo +=-.75*c,"Up",i,"Up",j
			# up down
			ampo +=-.25*c,"Up",i,"Down",j
			# down up 
			ampo +=-.25*c,"Down",i,"Up",j
			# down down
			ampo +=-.75*c,"Down",i,"Down",j

		end
		#Electric field#
		if !iszero(Estrength)
			ampo += -cos(angle)*Estrength,"X",i
			ampo += -sin(angle)*Estrength*fac2,"Y",i
		end
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
    #close(f)
    obs = DemoObserver(e_cutoff)
    energy,psi = dmrg(H,psi0,sweeps,observer=obs, outputlevel=0)
	mps_out=h5open(string("dmrg_logs/", dmrg_output_filename,string(round(g,digits=3))),"w")
	write(mps_out,"MPS",psi)
	close(mps_out)

	psi0 = copy(psi)

	f=open("dmrg_logs/log","a")
	println(f)
	maxbond=maxlinkdim(psi)
	println(f,"Max. bond dimension: ",maxbond)
	println(f)
	@printf(f,"Final energy = %.8f \n",energy)
	println(f)
	#close(f)

	if Nstates != 0 #&& (g==1.1 || g==0)
	#if Nstates != 0 
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

			nsweeps = 200
			maxdim = [200]
			cutoff = [1E-10]
			#noise = [1E-6]
			weight=30.

			psi_init = randomMPS(sites,Nbonds)
			#energy ,psi = dmrg(H,wavefunction[1:istates],initial_states[istates] ,sweeps, observer=obs, outputlevel=0)
			energy ,psi = dmrg(H,wavefunction[1:istates],psi_init ; nsweeps, maxdim, cutoff,weight)

			Hψ=apply(H, psi)
			normalize!(Hψ)
			println(inner(psi, Hψ))

			global initial_states[istates] = psi
			f2=open("dmrg_logs/log","a")
			maxbond=maxlinkdim(psi)
			println(f,"Max. bond dimension: ",maxbond)
			println(f2)
			println(f2,"Final energy "*string(istates)*". excited state= "*string(round(energy,digits=12))*"\n")
			println(f2)
			close(f2)
			wavefunction[istates+1]=psi
			append!(energies,energy)
			#Calculate ground state spectrum
			norm_mu=sqrt(inner(psi0, psi0)*inner(wavefunction[istates+1],wavefunction[istates+1]))
			intX=inner(wavefunction[istates+1],  Apply(muX, psi0))/norm_mu
			intY=inner(wavefunction[istates+1],  Apply(muY, psi0))/norm_mu
			println(f_int,g," ",energy," ",intX," ",intY,"\n")
		end

		# PN addition: orthoginalize excited states; even if that are not converged, they for a basis to diagonalize H
		Hmatrix=zeros(Float64,(Nstates+1, Nstates+1))
		Smatrix=zeros(Float64,(Nstates+1, Nstates+1))
		muXmatrix=zeros(Float64,(Nstates+1, Nstates+1))
		spectrum1=zeros(Float64,Nstates+1)
		spectrum2=zeros(Float64,Nstates+1)

		for i=1:Nstates+1
			normi=inner(wavefunction[i],wavefunction[i])
			for j=1:Nstates+1
				normj=inner(wavefunction[j],wavefunction[j])
				Hmatrix[i,j]=real(inner(wavefunction[i],  Apply(H, wavefunction[j]))/normi/normj)
				muXmatrix[i,j]=real(inner(wavefunction[i],  Apply(muX, wavefunction[j]))/normi/normj)
				Smatrix[i,j]=real(inner(wavefunction[i],wavefunction[j])/normi/normj)
			end
		end
		F = eigen(Hmatrix, Smatrix)

		spectrum1x=zeros(Float64,Nstates+1)
		spectrum1y=zeros(Float64,Nstates+1)
		spectrum2x=zeros(Float64,Nstates+1)
		spectrum2y=zeros(Float64,Nstates+1)

		Sinv=real(inv(Smatrix))
		Sinvsqrt=real(sqrt(Sinv))



		Hprime=real(copy(Hmatrix))
		Hp=real(copy(Hmatrix))
		mul!(Hprime,Hmatrix,Sinvsqrt)
		mul!(Hp,Sinvsqrt,Hprime)
		muprime=real(copy(muXmatrix))
		muXmatrixS=real(copy(muXmatrix))
		mul!(muprime,muXmatrix,Sinvsqrt)
		mul!(muXmatrixS,Sinvsqrt,muprime)
	
		Fp = eigen(Hp)
		Fpr=real(Fp.vectors)

		#display(heatmap(real(Smatrix)))
		#avefig("fig.png")

		sum1x=0.
		sum2x=0.
		for i=1:Nstates+1
			norm_mu1=real(sqrt(dot(Fp.vectors[:,1],Fp.vectors[:,1])*dot(Fp.vectors[:,i],Fp.vectors[:,i])))
			norm_mu2=real(sqrt(inner(wavefunction[1], wavefunction[1])*inner(wavefunction[i],wavefunction[i])))
			spectrum1x[i]=dot(Fpr[:,i], real(muXmatrixS), Fpr[:,1])/norm_mu1
			spectrum2x[i]=inner(wavefunction[i],  Apply(muX, wavefunction[1]))/norm_mu2
			sum1x+=spectrum1x[i]^2
			sum2x+=spectrum2x[i]^2
		end
		println(sum1x,' ',sum2x)

		text_energy=" "
		text_ent_vN=" "
		text_ent_R=" "
		text_xcorr=" "
		text_ycorr=" "
		text_corr=" "
		text_mux=" "
		text_muy=" "
		for istates=1:Nstates+1
			#FS.values[istates]=1.0/FS.values[istates]
			println(Fp.values[istates],' ',spectrum1x[istates]^2,' ',energies[istates],' ',spectrum2x[istates]^2,' ',Hmatrix[istates,istates],' ',Smatrix[1,istates])

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

		write_text("dmrg_logs/energy.txt",g,text_energy)
		write_text("dmrg_logs/entropy.txt",g,text_ent_vN)
		write_text("dmrg_logs/xcorr.txt",g,text_xcorr)
		write_text("dmrg_logs/ycorr.txt",g,text_ycorr)
		write_text("dmrg_logs/corr.txt",g,text_corr)
		write_text("dmrg_logs/mux.txt",g,text_mux)
		write_text("dmrg_logs/muy.txt",g,text_muy)

	else
		#Write energy to file#
		f=open("dmrg_logs/energy.txt","a")
		println(f,round(g,digits=4)," ",energy)
		#close(f)
	
		#Write entropy and Schmidt values to file#
		SvN,Renyi,Svalues = vN_entropy(psi,mbond)
	
		write_output("dmrg_logs/entropy_vN.txt",g,SvN)
		write_output("dmrg_logs/entropy_Renyi.txt",g,Renyi)
		write_output("dmrg_logs/schmidt_values.txt",g,Svalues)
		
		#Write nearest-neighbour correlation to file#
		xcorr,ycorr = correlation(psi,Nsites,Nspec,evod,X,Y)
		
		write_output("dmrg_logs/xcorr.txt",g,xcorr)
		write_output("dmrg_logs/ycorr.txt",g,ycorr)
		write_output("dmrg_logs/corr.txt",g,(xcorr+ycorr)/(Nsites-1))
		
		
		#Write polarization to file#
		MuX,MuY = polarization(psi,Nsites,Nspec,evod,X,Y)
	
		write_output("dmrg_logs/mux.txt",g,MuX)
		write_output("dmrg_logs/muy.txt",g,MuY)
	end

end
end
f=open("dmrg_logs/log","a")
println(f)
println(f,"Calculation finished.")
#close(f)
close(f_int)
