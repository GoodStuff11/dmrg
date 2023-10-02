module observables
using LinearAlgebra
import KrylovKit

export correlation, vN_entropy

function correlation(state;Nsites, mmax)
    ## 1/(N-1) sum_1^(N-1) 1/2(U_i D_(i+1) + D_i U_(i+1))
    dim = 2*mmax + 1
    final_state = zeros(ComplexF64, size(state))
    for i in 1:Nsites-1
        j = i+1
        for k in 1:dim^Nsites
            mk_i = ((k-1)%dim^i) ÷ dim^(i-1)
            mk_j = ((k-1)%dim^j) ÷ dim^(j-1)

            # up down
            if (mk_i < dim - 1) && (mk_j > 0)
                final_state[k + dim^(i-1) - dim^(j-1)] += 0.5*state[k]
            end

            # down up
            if (mk_i > 0) && (mk_j < dim - 1)
                final_state[k - dim^(i-1) + dim^(j-1)] += 0.5*state[k]
            end
        end
    end
    return real(dot(conj.(state), final_state)/(Nsites - 1))
end

function vN_entropy(state; Nsites, mmax, split)
    dim = 2*mmax + 1
    state_matrix =  reshape(state,(dim^(Nsites-split),dim^split))'
    vector_size = dim^split # dim^(Nsites-split)
    n_singular_values = min(dim^(Nsites-split), dim^split)
    vals, _, _, info = KrylovKit.svdsolve(state_matrix, ones((vector_size,)), min(n_singular_values,30), :LR)

    SvN = 0.0
    renyi = 0.0
    p = vals.^2
    p ./= sum(p)
    for n in eachindex(vals)
        SvN -= p[n] * log(p[n])
        renyi += p[n]^2
    end
    renyi = - 0.5*log(renyi)
    # println(renyi)
    return SvN, renyi
end
end