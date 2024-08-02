module observables
using LinearAlgebra
import KrylovKit
using SparseArrays

export correlation, vN_entropy, reflection_projection, inversion_projection_m, parity_projection_m, inversion_projection_dvr

function correlation(state::Vector{T};Nsites, dim, axis="X") where {T<:Number}
    ## 1/(N-1) sum_1^(N-1) 1/2(U_i D_(i+1) + D_i U_(i+1))
    c = if axis=="X" 1 else -1 end
    final_state = zeros(T, size(state))
    for i in 1:Nsites-1
        j = i+1
        for k in 1:dim^Nsites
            mk_i = ((k-1)%dim^i) ÷ dim^(i-1)
            mk_j = ((k-1)%dim^j) ÷ dim^(j-1)

            # up down
            if (mk_i < dim - 1) && (mk_j > 0)
                final_state[k + dim^(i-1) - dim^(j-1)] += state[k]
            end

            # down up
            if (mk_i > 0) && (mk_j < dim - 1)
                final_state[k - dim^(i-1) + dim^(j-1)] += state[k]
            end

            # down down
            if (mk_i > 0) && (mk_j > 0) 
                final_state[k - dim^(i-1) - dim^(j-1)] += c*state[k]
            end

            # up up
            if (mk_i < dim - 1) && (mk_j < dim - 1)
                final_state[k + dim^(i-1) + dim^(j-1)] += c*state[k]
            end
        end
    end
    return real(dot(state, final_state)/4)
end

function vN_entropy(state; Nsites, dim, split)
    state_matrix =  sparse(reshape(state,(dim^(Nsites-split),dim^split))')
    vector_size = dim^split # dim^(Nsites-split)
    p = abs.(KrylovKit.eigsolve(state_matrix*state_matrix', vector_size, 10, :LR)[1])
    # n_singular_values = min(dim^(Nsites-split), dim^split)
    # using a random initial vector because that is unlikely to cause svdsolve to crash
    # vals, _, _, _ = KrylovKit.svdsolve(state_matrix, vector_size, min(n_singular_values,30), :LR)

    SvN = 0.0
    renyi = 0.0
    # p = vals.^2
    # p ./= sum(p)
    for n in eachindex(p)
        SvN -= p[n] * log(p[n])
        renyi += p[n]^2
    end
    renyi = - 0.5*log(renyi)
    # println(renyi)
    return SvN, renyi, p
end

function projection_helper(state, projected_state, index, transformed_index,parity)
    if parity == "even"
        if transformed_index == index
            # |abcba><abcba|
            projected_state[index] += state[index]
        else
            # (|abcde><abcde| + |f(abcde)><abcde|)/2
            projected_state[index] += state[index]/2
            projected_state[transformed_index] += state[index]/2
        end
    else
        # (|abcde><abcde| - |f(abcde)><abcde|)/2
        projected_state[index] += state[index]/2
        projected_state[transformed_index] -= state[index]/2
    end
end

function reflection_projection(state; Nsites, dim, parity="even")
    projected_state = zero(state)

    for index in 0:length(state)-1
        reflected_index = sum([((index÷dim^k)%dim)*dim^(Nsites-k-1) for k in 0:Nsites-1])
        # println(index," ",reflected_index, " ", dim^Nsites)
        projection_helper(state, projected_state, index+1, reflected_index+1, parity)
    end

    return projected_state
end


function inversion_projection_m(state; Nsites, dim, parity="even")
    # all m -> all -m
    @assert dim % 2 == 1
    projected_state = zero(state)

    for index in 0:length(state)-1
        inverted_index = sum([(dim-(index÷dim^k)%dim-1)*dim^k for k in 0:Nsites-1])
        projection_helper(state, projected_state, index+1, inverted_index+1, parity)
    end

    return projected_state
end

function parity_projection_m(state; Nsites, dim, parity="even")
    # sum m = even/odd
    mmax = dim ÷ 2
    is_even = (dim+1)%2
    projected_state = zero(state)

    parity_increment = 0
    if parity == "even"
        parity_increment = 1
    end
    index_zero_odd = ((mmax + is_even)*Nsites)%2
    for index in eachindex(state)
        if (index + parity_increment) % 2 == index_zero_odd
            projected_state[index] = state[index]
        end
    end

    return projected_state
end

function inversion_projection_dvr(state; Nsites, dim, parity="even")
    # updown: phi -> -phi

    projected_state = zero(state)

    for index in 0:length(state)-1
        inverted_index = sum([((dim - (index÷dim^k)%dim)%dim)*dim^k for k in 0:Nsites-1])
        projection_helper(state, projected_state, index+1, inverted_index+1, parity)
    end

    return projected_state
end

end   

