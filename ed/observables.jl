module observables
using LinearAlgebra
import KrylovKit

export correlation, vN_entropy, even_reflection_projection, odd_reflection_projection, even_inversion_projection, parity_projection

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
    # using a random initial vector because that is unlikely to cause svdsolve to crash
    vals, _, _, info = KrylovKit.svdsolve(state_matrix, rand((vector_size,)), min(n_singular_values,30), :LR)

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

function even_reflection_projection(state; Nsites, mmax)
    dim = 2*mmax + 1
    projected_state = zero(state)

    for index in 0:length(state)-1
        reflected_index = sum([((index÷dim^k)%dim)*dim^(Nsites-k-1) for k in 0:Nsites-1])
        # println(index," ",reflected_index, " ", dim^Nsites)
        if reflected_index == index
            # |abcba><abcba|
            projected_state[index+1] = state[index+1]
        else
            # (|abcde><abcde| + |edcba><abcde|)/2
            projected_state[index+1] += state[index+1]/2
            projected_state[reflected_index+1] += state[index+1]/2
        end
    end

    return projected_state
end

function odd_reflection_projection(state; Nsites, mmax)
    dim = 2*mmax + 1
    projected_state = zero(state)

    for index in 0:length(state)-1
        reflected_index = sum([((index÷dim^k)%dim)*dim^(Nsites-k-1) for k in 0:Nsites-1])
        # println(index," ",reflected_index, " ", dim^Nsites)
        if reflected_index != index
            # (|abcde><abcde| - |edcba><abcde|)/2
            projected_state[index+1] += state[index+1]/2
            projected_state[reflected_index+1] -= state[index+1]/2
        end
    end

    return projected_state
end

function even_inversion_projection(state; Nsites, mmax)
    # all m -> all -m
    dim = 2*mmax + 1
    projected_state = zero(state)

    for index in 0:length(state)-1
        inverted_index = sum([(dim-(index÷dim^k)%dim-1)*dim^k for k in 0:Nsites-1])
        # println(index," ",inverted_index, " ", dim^Nsites)
        if inverted_index == index
            # |abcba><abcba|
            projected_state[index+1] = state[index+1]
        else
            # (|abcde><abcde| + |edcba><abcde|)/2
            projected_state[index+1] += state[index+1]/2
            projected_state[inverted_index+1] += state[index+1]/2
        end
    end

    return projected_state
end

function parity_projection(state; Nsites, mmax, parity="even")
    # sum m = even/odd
    projected_state = zero(state)

    parity_increment = 0
    if parity == "even"
        parity_increment = 1
    end
    index_zero_odd = (mmax*Nsites)%2
    for index in 1:length(state)
        if (index + parity_increment) % 2 == index_zero_odd
            projected_state[index] = state[index]
        end
    end

    return projected_state
end


end   

