using ITensors
#using ITensorTDVP
using Printf
using Observers
using Printf
using HDF5
using StatsPlots
using CSV
using DataFrames
push!(LOAD_PATH,joinpath(pwd(),"dmrg_tdvp"))

using utility_funcs
using operators
using states
using dvr

import YAML


function main(data::Dict{Any, Any})
    mmax_min = data["mmax_min"]
    mmax_max = data["mmax_max"]
    mmax_diff = data["mmax_diff"]

    Nsites_min = data["Nsites_min"]
    Nsites_max = data["Nsites_max"]
    Nsites_diff = data["Nsites_diff"]
    Nbonds = data["Nbond"]
    Nsweep = data["Nsweep"]
    e_cutoff = data["ecutoff"]
    SVD_error = data["SVD"]
    gstart = data["gstart"]
    delta_g = data["dg"]
    Ng = data["Ng"]
    s_split = data["s_split"]
    pairs = data["pairs"]
    evod = data["states"]
    angle = data["angle"]*pi/180.0
    Estrength = data["strength"]
    Nstates = data["Nstates"]
    maxdim = data["maxdim"]
    output_path = data["outputpath"]

    set_logger(output_path)

    if evod == "dvr"
        fac1 = 1.0
        fac2 = 1.0
    else 
        fac1 = -1.0
        fac2 = 1.0im
    end

    #Define basis#
    if evod == "all"
        @info "all m-states are considered"
    elseif evod == "dvr"
        @info  "DVR-basis is used"
    else
        @info ("only "+evod+" m-states are considered")
    end


    if pairs == "nearest"
        @info "only nearest-neighbour interactions"
    elseif pairs == "allpairs"
        @info "all interactions"
    end

    @info """\n#################################################################################
        ##################################
        ####Calculate free rotor chain####
        ##################################
        """


    #Define output files#
    if Nstates == 0
        df = DataFrame(dmrg_duration=Float64[], bondsize=Int[], memory=Int[], 
                mmax=Int[],Nsites=Int[],interaction=String[],basis=String[],
                E_angle=[], E_strength=[],
                g=[],energy=[], SvN=[], renyi=[],mux=[],muy=[], xcorr=[],ycorr=[],corr=[])
    else
        df = DataFrame(dmrg_duration=Float64[], bondsize=Int[],memory=Int[],
                mmax=Int[],Nsites=Int[],interaction=String[], basis=String[],
                E_angle=[], E_strength=[],
                g=[],istate=Int[],energy=[], SvN=[], renyi=[],mux=[],muy=[], xcorr=[],ycorr=[],corr=[])
    end
    
    for mmax = mmax_min:mmax_diff:mmax_max
        for Nsites = Nsites_min:Nsites_diff:Nsites_max
            for i=1:Nstates+1
                create_file(joinpath(output_path,"schmidt_values_"*string(i)*".txt"),evod,mmax,Nsites)
            end
            @info @sprintf """\n#######################################
            ###########Basis information###########
            #######################################
            mmax = %i
            Number of sites: %i
            Dimension of local Hilbert space:: %i
            """ mmax Nsites 2*mmax+1
        
            @info @sprintf "Dimension of local Hilbert space for chosen m-states: %i" 2*mmax+1
            #Determine number of interaction pairs per starting site#
            Nsecond = zeros(Int64,(Nsites-1))
            for i=1:Nsites-1
                if pairs == "nearest"
                    Nsecond[i]=i+1
                elseif pairs == "allpairs"
                    Nsecond[i]=Nsites
                end
            end

            #Non-interacting rotors as initial guess#
            initial_states = nothing
            sites=siteinds("PlaRotor",Nsites; conserve_qns=false, dim=mmax)
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
            psi0 = randomMPS(sites,10) # global
            
            maxdim!(sweeps,10) # gradually increase states kept
            cutoff!(sweeps,SVD_error) # desired truncation error
            
            #Perform DMRG runs#
            obs = DemoObserver(e_cutoff)
            duration = @elapsed begin
                energy,psi = dmrg(H0,psi0,sweeps,observer=obs, outputlevel=0)
            end
            psi0 = psi #global

            maxbond=maxlinkdim(psi)
            # push!(bonddim_log, maxbond)
            @info @sprintf """\nMax. bond dimension: %i
            DMRG Duration: %f
            Final energy = %.8f \n
            Initial state calculated
            ###############################################################################
            """ maxbond duration energy

            if s_split == -1
                mbond = Nsites÷2
            else
                mbond = s_split # customize bond
            end

            for ig = 0:Ng-1
                duration_log = Vector{Float64}()
                bonddim_log = Vector{Int}()
                memory_log = Vector{Int}()

                g = gstart + ig*delta_g
                
                @info @sprintf("""\n##################################
                    ########g= %.3f ########
                    ##################################

                    ####DMRG calculation####
                    Construct MPO""", g
                )
                
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
                @info "Start DMRG run"

                obs = DemoObserver(e_cutoff)
                duration = @elapsed begin
                    memory = @allocated(tmp = dmrg(H,psi0,sweeps,observer=obs, outputlevel=0))
                    energy,psi = tmp
                    push!(memory_log, memory)
                end
                push!(duration_log, duration)
                mps_out=h5open(joinpath(output_path,"psi0_g"*string(round(g,digits=3))),"w")
                write(mps_out,"MPS",psi)
                close(mps_out)

                psi0 = copy(psi)

                maxbond=maxlinkdim(psi)
                push!(bonddim_log, maxbond)
                @info @sprintf """\nMax. bond dimension: %i
                DMRG Duration: %f
                Final energy = %.8f \n
                """ maxbond duration energy
                
                
                if Nstates > 0 #&& (g==1.1 || g==0)
                #if Nstates != 0 

                    energies = []
                    append!(energies,energy)
                    wavefunction = [psi for ii=1:Nstates+1]
                    if ig == 0
                        initial_states = [psi0 for ii=1:Nstates] #global
                    end

                    for istates=1:Nstates
                        if ig == 0
                            initial_states[istates] = randomMPS(sites,Nbonds) #global
                        else
                            maxdim!(sweeps,maxlinkdim(initial_states[istates]),Nbonds) # gradually increase states kept
                        end
                        cutoff!(sweeps,SVD_error) # desired truncation error

                        cutoff = [1E-10]
                        #noise = [1E-6]
                        weight=30.

                        psi_init = randomMPS(sites,Nbonds)
                        #energy ,psi = dmrg(H,wavefunction[1:istates],initial_states[istates] ,sweeps, observer=obs, outputlevel=0)
            
                        duration = @elapsed begin
                            memory = @allocated(tmp = dmrg(H,wavefunction[1:istates],psi_init; nsweeps=Nsweep, maxdim=maxdim, cutoff=cutoff,weight=weight))
                            push!(memory_log, memory)
                            energy ,psi = tmp
                        end
                        push!(duration_log, duration)
                        initial_states[istates] = psi # global

                        maxbond=maxlinkdim(psi)
                        push!(bonddim_log, maxbond)
                        @info @sprintf("""\nMax. bond dimension: %i
                            Final energy %i - excited state=%.12f\n
                            """, maxbond, istates, round(energy,digits=12))
                        
                        wavefunction[istates+1]=psi
                        append!(energies,energy)
                        #Calculate ground state spectrum
                        norm_mu=sqrt(inner(psi0, psi0)*inner(wavefunction[istates+1],wavefunction[istates+1]))
                        intX=inner(wavefunction[istates+1],  Apply(muX, psi0))/norm_mu
                        intY=inner(wavefunction[istates+1],  Apply(muY, psi0))/norm_mu
                        @info @sprintf("g=%f energy=%f intX=%f intY=%f", g,energy,intX,intY)
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
                    @info @sprintf("%f %f", sum1x, sum2x)


                    for istates=1:Nstates+1
                        #FS.values[istates]=1.0/FS.values[istates]
                        @info @sprintf("%f %f %f %f %f %f", Fp.values[istates],spectrum1x[istates]^2,energies[istates],spectrum2x[istates]^2,Hmatrix[istates,istates],Smatrix[1,istates])

                        #Calculate von-Neumann entropy and Schmidt coefficients#
                        SvN,Renyi,Svalues = vN_entropy(wavefunction[istates],mbond)
                    
                        write_output("schmidt_values_"*string(istates)*".txt",g,Svalues)
                        
                        #Calculate dipole correlations#
                        xcorr,ycorr = correlation(wavefunction[istates],Nsites,evod)
                        
                        #Calculate summed dipole moment and fluctuation#
                        MuX,MuY = polarization(wavefunction[istates],Nsites,evod)   
                        push!(df, [duration_log[istates], bonddim_log[istates], memory_log[istates],
                            mmax, Nsites, pairs, evod, 
                            angle, Estrength,
                            g, istates, energies[istates], SvN, Renyi, MuX, MuY, xcorr, ycorr, (xcorr+ycorr)/(Nsites-1)])
                    end    
                    
                
                else    
                    #Write entropy and Schmidt values to file#
                    SvN,Renyi,Svalues = vN_entropy(psi,mbond)
                    #Write nearest-neighbour correlation to file#
                    xcorr,ycorr = correlation(psi,Nsites,evod)
                    #Write polarization to file#
                    MuX,MuY = polarization(psi,Nsites,evod)

                    write_output(joinpath(output_path, "schmidt_values.txt"),g,Svalues)
                    
                    push!(df, [
                        duration_log[1], bonddim_log[1], memory_log[1],
                        mmax, Nsites, pairs, evod, 
                        angle, Estrength,
                        g, energy, SvN, Renyi, MuX, MuY, xcorr, ycorr, (xcorr+ycorr)/(Nsites-1)])
                end

            end
        end
    end

    CSV.write(joinpath(output_path, "data.csv"), df)
    @info "Calculation finished."
end

data = YAML.load_file("input_quick.yml")

if length(ARGS) > 0
    data["outputpath"] = joinpath(data["outputpath"],ARGS[1])
    if !isdir(data["outputpath"])
        mkdir(data["outputpath"])
    end
    # data["gstart"] = parse(Float64, ARGS[2])
    # data["dg"] = parse(Float64, ARGS[3])
    # data["Ng"] = parse(Int, ARGS[4])
    # data["mmax_min"] = parse(Int,ARGS[5])
    # data["Nsites"] = parse(Int,ARGS[6])
    # data["mbond"] = data["Nsites"] ÷ 2
end

main(data)
