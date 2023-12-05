import KrylovKit
using CSV
using DataFrames
using LinearAlgebra

push!(LOAD_PATH,joinpath(pwd(),"ed"))
using observables

function Hamiltonian(state; Nsites, mmax, g, angle, Estrength, pairs)
    dim = 2*mmax + 1
    final_state = zeros(ComplexF64, size(state))
    for i in 1:Nsites
        # Ti
        for k in 1:dim^Nsites # iterate over all states
            mk = ((k-1)%dim^i) ÷ dim^(i-1) # m value at atom from 0 to dim-1
            m2 = (mk - (dim-1)/2)^2
            
            # diagonal
            final_state[k] += state[k]*m2

            # upper diagonal
            if mk < dim - 1
                final_state[k+dim^(i-1)] += 0.5*(
                    -cos(angle)*Estrength - sin(angle)*Estrength*1im
                    )*state[k]
            end
            # lower diagonal
            if mk > 0
                final_state[k-dim^(i-1)] += 0.5*(
                    -cos(angle)*Estrength + sin(angle)*Estrength*1im
                    )*state[k]
            end

        end
        # Vij
        if g == 0
            continue
        end

        site_pairing = (pairs == "nearest") ? min(i+1, Nsites) : Nsites
        for j in (i+1):site_pairing
            c = g/abs(j-i)^3
            for k in 1:dim^Nsites
                mk_i = ((k-1)%dim^i) ÷ dim^(i-1)
                mk_j = ((k-1)%dim^j) ÷ dim^(j-1)
                # up up
                if (mk_i < dim - 1) && (mk_j < dim - 1)
                    final_state[k + dim^(i-1) + dim^(j-1)] -= 0.75*c*state[k]
                end
                # up down
                if (mk_i < dim - 1) && (mk_j > 0)
                    final_state[k + dim^(i-1) - dim^(j-1)] -= 0.25*c*state[k]
                end

                # down up
                if (mk_i > 0) && (mk_j < dim - 1)
                    final_state[k - dim^(i-1) + dim^(j-1)] -= 0.25*c*state[k]
                end
                # down down
                if (mk_i > 0) && (mk_j > 0)
                    final_state[k - dim^(i-1) - dim^(j-1)] -= 0.75*c*state[k]
                end
            end
        end
    end
    return final_state
end


function main()
    g = 1.2
    pairs = "nearest"
    Nsites = 2
    @show Nsites
    mmax = 1
    initial_state = rand(ComplexF64, (2*mmax+1)^Nsites)

    state1 = odd_reflection_projection(even_inversion_projection(initial_state; 
            Nsites=Nsites, mmax=mmax);
            Nsites=Nsites, mmax=mmax)
    state2 = even_inversion_projection(odd_reflection_projection(initial_state; 
            Nsites=Nsites, mmax=mmax);
            Nsites=Nsites, mmax=mmax)
    println(dot(state1,state2)/(sqrt(dot(state1,state1))*sqrt(dot(state2,state2))))
    # for mmax in 2
    #     println("mmax",mmax)
    #     Nsites = 6
    #     println("Nsites",Nsites)
    #     prev_vecs = nothing
    #     for g in 0.01:0.05:3
    #         energies = 20
    #         Estrength = 0
    #         angle = 90
    #         outputpath = "/home/jkambulo/projects/def-pnroy/jkambulo/dmrg/output_data"
            
    #         filename = joinpath(outputpath, "deleteme.csv")
    #         if isfile(filename)
    #             df = DataFrame(CSV.File(filename))
    #         else
    #             columns = Dict("t_diagonalization"=>Float64[],"t_metrics"=>Float64[],"pairs"=>String[], "memory"=>Int[], "mmax"=>Int[], 
    #                             "Nsites"=>Int[], "g"=>Float64[],"Estrength"=>Float64[], "angle"=>Float64[], )
    #             for i in 1:energies
    #                 columns["E$i"] = Float64[]
    #                 columns["correlation$i"] = Float64[]
    #                 columns["SvN$i"] = Float64[]
    #                 columns["reflection_symmetry$i"] = Float64[]
    #                 columns["inversion_symmetry$i"] = Float64[]
    #                 columns["parity$i"] = Float64[]
    #                 columns["prev_corresponding_energy$i"] = Int[]

    #             end
    #             df = DataFrame(columns)
    #         end
    #         @elapsed begin end
    #         t1 = @elapsed begin
    #             initial_state = ones(ComplexF64, (2*mmax+1)^Nsites)
    #             H(x) = Hamiltonian(x; Nsites=Nsites, mmax=mmax, g=g, angle=angle, Estrength=Estrength, pairs=pairs)
    #             memory = @allocated(tmp = KrylovKit.eigsolve(H, initial_state, energies, :SR))
    #             vals, vecs, info = tmp
    #         end
    #         data = Dict("memory"=>memory, "mmax"=>mmax, "Nsites"=>Nsites,
    #             "g"=>g,"Estrength"=>Estrength, "angle"=>angle, "pairs"=>pairs)
    #         t2 = @elapsed begin
    #             for i in 1:energies
    #                 data["E$i"] = real(vals[i])
    #                 data["correlation$i"] = correlation(vecs[i]; Nsites=Nsites, mmax=mmax)
    #                 SvN, _ = vN_entropy(vecs[i]; Nsites=Nsites, mmax=mmax, split=Nsites ÷ 2)
    #                 data["SvN$i"] = SvN                
    #                 data["reflection_symmetry$i"]= dot(vecs[i], even_reflection_projection(vecs[i]; Nsites=Nsites, mmax=mmax))
    #                 data["inversion_symmetry$i"]= dot(vecs[i], even_inversion_projection(vecs[i]; Nsites=Nsites, mmax=mmax))
    #                 data["parity$i"] = dot(vecs[i], parity_projection(vecs[i]; Nsites=Nsites, mmax=mmax, parity="even"))
    #                 if prev_vecs == nothing
    #                     data["prev_corresponding_energy$i"] = i
    #                 else
    #                     diffs = [abs(dot(vecs[i], v1)) for v1 in prev_vecs]
    #                     data["prev_corresponding_energy$i"] = argmax(diffs)
    #                     if diffs[data["prev_corresponding_energy$i"]] <0.5
    #                         data["prev_corresponding_energy$i"] = -1
    #                     end
    #                 end
    #                 # if prev_vecs != nothing
    #                 #     println(diffs," ", data["prev_corresponding_energy$i"])
    #                 # end
    #             end
    #         end
    #         data["t_diagonalization"] = t1
    #         data["t_metrics"] = t2
    #         # println(names(df))
    #         # println(data)
    #         prev_vecs = vecs[1:energies]
    #         push!(df, data)
    #         CSV.write(filename, df)
    #     end
    # end
end

main()