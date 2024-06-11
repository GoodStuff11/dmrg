using Symbolics
using LinearAlgebra
using IterTools


function dipole_chain_hamiltonian(Nsites, g)
    @variables x[1:Nsites]
    hamiltonian = g*sum(sin(ϕ_1)*sin(ϕ_2) - 2*cos(ϕ_1)*cos(ϕ_2) for (ϕ_1, ϕ_2) in partition(x,2,1))
    initial_conditions = Dict(y => 0 for y in x)
    return x, hamiltonian, initial_conditions
end



function get_first_order_eigenvalues(args; n_eigenvalues=10)
    x, hamiltonian, initial_conditions = args
    ground_energy = substitute(hamiltonian, initial_conditions).val
    a = substitute(Symbolics.hessian(hamiltonian,[x...]), initial_conditions)
    a = map(x-> x.val, a)

    λ, _ = eigen(a)

    eigenvalues = []
    for n in Tuple.(CartesianIndices(Tuple(0:n_eigenvalues for i=1:Nsites)))
        push!(eigenvalues,sum((n .+ 1/2) .* sqrt.(2*λ)) + ground_energy)
    end
    return sort(eigenvalues)[1:n_eigenvalues]
end
# get_first_order_eigenvalues(dipole_chain_hamiltonian(5, 1.1))