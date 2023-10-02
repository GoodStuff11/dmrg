using ITensors
#using ITensorTDVP
using Observers
using Printf
using HDF5
using StatsPlots
using CSV
using DataFrames
push!(LOAD_PATH,joinpath(pwd(),"dmrg_tdvp"))

using matrices
using states
using expectations
using dvr
using utility_funcs

import YAML


function main(data::Dict{Any, Any})
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

    logger = Logger(output_path, "log")
    logger_int = Logger(output_path, "int.txt")
    logger < [
        """
        #######################################
        ###########Basis information###########
        #######################################""",
        ("mmax= ",mmax),
        ("Number of sites: ",Nsites),
        ("Dimension of local Hilbert space: ",2*mmax+1)
    ]

    if evod == "dvr"
        fac1 = 1.0
        fac2 = 1.0
    else 
        fac1 = -1.0
        fac2 = 1.0im
    end

    #Calculate kinetic matrix and x operator#
    Ttmp = kinetic(mmax)
    Xtmp = Xoperator(mmax)
    Ytmp = Yoperator(mmax)
    Uptmp = Upoperator(mmax)
    Downtmp = Downoperator(mmax)
    #Define basis#
    if evod == "all"
        global T = Ttmp
        global X = Xtmp
        global Y = Ytmp
        global Up = Uptmp
        global Down = Downtmp

        #hilbert space size
        Nspec=size(T,1)
        logger < "all m-states are considered"
    elseif evod == "dvr"
        tmp1,tmp2,tmp3 = exp_dvr(mmax)
        global T = tmp1
        global X = tmp2
        global Y = tmp3
        Nspec=size(T,1)
        logger < "DVR-basis is used"
    else
        global T = symmetry(Ttmp,evod,mmax)
        global X = symmetry(Xtmp,evod,mmax)
        global Y = symmetry(Ytmp,evod,mmax)
        logger < ("only "+evod+" m-states are considered")
        Nspec=size(T,1)
    end
    logger < ("Dimension of local Hilbert space for chosen m-states: ",Nspec)
    if pairs == "nearest"
        logger < "only nearest-neighbour interactions"
    elseif pairs == "allpairs"
        logger < "all interactions"
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

    logger <
        """\n
        #################################################################################
        ##################################
        ####Calculate free rotor chain####
        ##################################
        """


    #Define output files#
    if Nstates == 0
        df = DataFrame(mmax=Int[],Nsites=Int[],interaction=String[],basis=String[],
                E_angle=[], E_strength=[],
                g=[],energy=[], vN=[], renyi=[],mux=[],muy=[], xcorr=[],ycorr=[],corr=[])
    else
        df = DataFrame(mmax=Int[],Nsites=Int[],interaction=String[], basis=String[],
                E_angle=[], E_strength=[],
                g=[],istate=Int[],energy=[], vN=[], renyi=[],mux=[],muy=[], xcorr=[],ycorr=[],corr=[])
    end


    for i=1:Nstates+1
        create_file(joinpath(output_path,"schmidt_values_"*string(i)*".txt"),evod,mmax,Nsites)
    end

    let
    include("./dmrg_tdvp/operators.jl")
    include("./dmrg_tdvp/observer.jl")
    for ig = 0:Ng-1
        #Non-interacting rotors as initial guess#
        if ig == 0
            sites=siteinds("PlaRotor",Nsites; conserve_qns=false)
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
            duration = @elapsed begin
                energy,psi = dmrg(H0,psi0,sweeps,observer=obs, outputlevel=0)
            end
            global psi0 = psi

            maxbond=maxlinkdim(psi)
            logger < [("Max. bond dimension: ",maxbond),
                    ("DMRG Duration: ",duration),
                    @sprintf("Final energy = %.8f \n\n",energy),
                    "Initial state calculated",
                    "###############################################################################"]
        end    
        g = gstart + ig*delta_g
        
        logger < [
            "##################################",
            ("########g= ",g," ########"),
            """##################################

            ####DMRG calculation####
            Construct MPO"""
        ]
        

        # setup MPO for the correct Hamiltonian
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
        logger < "Start DMRG run"

        obs = DemoObserver(e_cutoff)
        duration = @elapsed begin
            energy,psi = dmrg(H,psi0,sweeps,observer=obs, outputlevel=0)
        end
        mps_out=h5open(joinpath(output_path,"psi0_g"*string(round(g,digits=3))),"w")
        write(mps_out,"MPS",psi)
        close(mps_out)

        psi0 = copy(psi)

        maxbond=maxlinkdim(psi)
        logger < ["",
            ("Max. bond dimension: ",maxbond),
            ("DMRG Duration: ",duration),
            @sprintf("Final energy = %.8f \n", energy),
            ""
            ]
        
        if Nstates > 0 #&& (g==1.1 || g==0)
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
                global initial_states[istates] = psi

                maxbond=maxlinkdim(psi)
                logger < [
                    ("Max. bond dimension: ",maxbond),"",
                    "Final energy "*string(istates)*". excited state= "*string(round(energy,digits=12))*"\n",
                    ""
                    ]
                wavefunction[istates+1]=psi
                append!(energies,energy)
                #Calculate ground state spectrum
                norm_mu=sqrt(inner(psi0, psi0)*inner(wavefunction[istates+1],wavefunction[istates+1]))
                intX=inner(wavefunction[istates+1],  Apply(muX, psi0))/norm_mu
                intY=inner(wavefunction[istates+1],  Apply(muY, psi0))/norm_mu
                logger_int < (g," ",energy," ",intX," ",intY,"\n")
            end
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
            logger < (sum1x,' ',sum2x)


            for istates=1:Nstates+1
                #FS.values[istates]=1.0/FS.values[istates]
                logger < (Fp.values[istates],' ',spectrum1x[istates]^2,' ',energies[istates],' ',spectrum2x[istates]^2,' ',Hmatrix[istates,istates],' ',Smatrix[1,istates])

                #Calculate von-Neumann entropy and Schmidt coefficients#
                SvN,Renyi,Svalues = vN_entropy(wavefunction[istates],mbond)
            
                write_output("schmidt_values_"*string(istates)*".txt",g,Svalues)
                
                #Calculate dipole correlations#
                xcorr,ycorr = correlation(wavefunction[istates],Nsites,Nspec,evod,X,Y)
                
                #Calculate summed dipole moment and fluctuation#
                MuX,MuY = polarization(wavefunction[istates],Nsites,Nspec,evod,X,Y)   
                push!(df, [
                    mmax, Nsites, pairs, evod, 
                    angle, Estrength,
                    g, istates, energies[istates], SvN, Renyi, MuX, MuY, xcorr, ycorr, (xcorr+ycorr)/(Nsites-1)])
            end    
            
        
        else    
            #Write entropy and Schmidt values to file#
            SvN,Renyi,Svalues = vN_entropy(psi,mbond)
            #Write nearest-neighbour correlation to file#
            xcorr,ycorr = correlation(psi,Nsites,Nspec,evod,X,Y)
            #Write polarization to file#
            MuX,MuY = polarization(psi,Nsites,Nspec,evod,X,Y)

            write_output(joinpath(output_path, "schmidt_values.txt"),g,Svalues)
            
            push!(df, [
                mmax, Nsites, pairs, evod, 
                angle, Estrength,
                g, energy, SvN, Renyi, MuX, MuY, xcorr, ycorr, (xcorr+ycorr)/(Nsites-1)])
        end

    end

    end

    CSV.write(joinpath(output_path, "data.csv"), df)
    logger < "Calculation finished."
end

data = YAML.load_file("input_quick.yml")

if length(ARGS) > 0
    data["outputpath"] = joinpath(data["outputpath"],args[1])
    data["gstart"] = parse(Float64, args[2])
    data["dg"] = parse(Float64, args[3])
    data["Ng"] = parse(Int, args[4])
    data["m"] = parse(Int,args[5])
end

main(data)