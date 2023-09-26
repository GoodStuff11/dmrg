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
    g = 30
    Nsites = 6
    for g in 0.1:0.1:30
        for mmax in 1:5
            energies = 20
            Estrength = 0
            angle = 90
            outputpath = "C:\\Users\\jonat\\OneDrive\\Documents\\programming\\AnacondaProjects\\PHYS437A\\dmrg\\output_data"
            filename = joinpath(outputpath, "ED_benchmark.csv")
            if isfile(filename)
                df = DataFrame(CSV.File(filename))
            else
                columns = Dict("t"=>Float64[], "memory"=>Int[], "mmax"=>Int[], 
                                "Nsites"=>Int[], "g"=>Float64[],"Estrength"=>Float64[], 
                                "angle"=>Float64[], "correlation"=>Float64[])
                for i in 1:energies
                    columns["E$i"] = Float64[]
                end
                df = DataFrame(columns)
            end
            @elapsed begin
            end
            t = @elapsed begin
                initial_state = ones(ComplexF64, (2*mmax+1)^Nsites)
                H(x) = Hamiltonian(x; Nsites=Nsites, mmax=mmax, g=g, angle=angle, Estrength=Estrength, pairs="nearest")
                memory = @allocated(tmp = KrylovKit.eigsolve(H, initial_state, energies, :SR; krylovdim = 30))
                vals, vecs, info = tmp

                ground_corr = correlation(vecs[1]; Nsites=Nsites, mmax=mmax)
            end
            data = Dict("t"=>t, "memory"=>memory, "mmax"=>mmax, "Nsites"=>Nsites,
                        "g"=>g,"Estrength"=>Estrength, "angle"=>angle, "correlation"=>ground_corr)
            for i in 1:energies
                data["E$i"] = real(vals[i])
            end
            push!(df, data)
            CSV.write(filename, df)

        end
    end
end

main()