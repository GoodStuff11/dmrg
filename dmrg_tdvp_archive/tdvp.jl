using ITensors
using ITensorTDVP
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

    if !isdir(data["outputpath"])
        mkdir(data["outputpath"])
    end

    dt = data["dt"]
    N_timesteps = data["N_timesteps"]

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
    df = DataFrame(bondsize=Int[],
            mmax=Int[],Nsites=Int[],interaction=String[], basis=String[],
            E_angle=[], E_strength=[], g=[],
            autocorrelation=[],energy=[], SvN=[], renyi=[],mux=[],muy=[], xcorr=[],ycorr=[],corr=[], dt=[])

    
    for mmax = mmax_min:mmax_diff:mmax_max
        for Nsites = Nsites_min:Nsites_diff:Nsites_max
            # for i=1:Nstates+1
            #     create_file(joinpath(output_path,"schmidt_values_"*string(i)*".txt"),evod,mmax,Nsites)
            # end
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
            sites=siteinds("PlaRotor",Nsites; conserve_qns=false, dim=2*mmax+1)
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
                end
                # mps_out=h5open(joinpath(output_path,"psi0_g"*string(round(g,digits=3))),"w")
                # write(mps_out,"MPS",psi)
                # close(mps_out)

                maxbond=maxlinkdim(psi)
                @info @sprintf """\nMax. bond dimension: %i
                DMRG Duration: %f
                Final energy = %.8f \n
                """ maxbond duration energy

                start_psi = apply(muX, psi)
                energy = inner(psi, apply(H, start_psi))
                SvN, Renyi, Svalues = vN_entropy(psi, mbond)
                xcorr,ycorr = correlation(psi,Nsites,evod)
                polarizationX,polarizationY = polarization(psi,Nsites,evod)
                push!(df, [
                    maxbond,
                    mmax, Nsites, pairs, evod, 
                    angle, Estrength, g, 
                    1, energy, SvN, Renyi, polarizationX, polarizationY, xcorr, ycorr, 
                    (xcorr+ycorr)/(Nsites-1), 0
                ])
                t_values = range(dt, N_timesteps*dt, step=dt)
                
                t = @elapsed begin
                    for i=1:N_timesteps
                        println(i)
                        ψ_td = tdvp(H, -t_values[i]*1im, start_psi);
                        overlap = inner(psi,ψ_td)
                        energy = inner(ψ_td, apply(H, ψ_td))
                        SvN, Renyi, Svalues = vN_entropy(ψ_td, mbond)
                        xcorr,ycorr = correlation(ψ_td,Nsites,evod)
                        polarizationX,polarizationY = polarization(ψ_td,Nsites,evod)
                        push!(df, [
                            maxbond,
                            mmax, Nsites, pairs, evod, 
                            angle, Estrength, g, 
                            overlap, energy, SvN, Renyi, polarizationX, polarizationY, xcorr, ycorr, 
                            (xcorr+ycorr)/(Nsites-1), t_values[i]
                        ])
                    end
                end
                println("TDVP time: ", t)
                # write_output(joinpath(output_path, "schmidt_values.txt"),g,Svalues)
                
            end
        end
    end

    CSV.write(joinpath(output_path, "data.csv"), df)
    @info "Calculation finished."
end

data = YAML.load_file("input_quick.yml")

if length(ARGS) > 0
    data["outputpath"] = joinpath(data["outputpath"],ARGS[1])
    # data["gstart"] = parse(Float64, ARGS[2])
    # data["dg"] = parse(Float64, ARGS[3])
    # data["Ng"] = parse(Int, ARGS[4])
    # data["mmax_min"] = parse(Int,ARGS[5])
    # data["Nsites"] = parse(Int,ARGS[6])
    # data["mbond"] = data["Nsites"] ÷ 2
end

main(data)
