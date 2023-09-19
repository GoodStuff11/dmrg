import KrylovKit
using ITensors

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
    @time begin
        Nsites = 10
        mmax = 1
        initial_state = ones(ComplexF64, (2*mmax+1)^Nsites)
        H(x) = Hamiltonian(x; Nsites=Nsites, mmax=mmax, g=0.1, angle=90, Estrength=0, pairs="nearest")
        vals, vecs, info = KrylovKit.eigsolve(H, initial_state, 30, :SR; krylovdim = 30)
        print(vals)
        print(info)
    end
end

main()