function generate_initial_state(sites; parity_symmetry_type="none", inversion_symmetry_type="none")
	# generates an initial random MPS with desired parity and inversion symmetries. 
    Nsites = length(sites)
    Nspec = dim(sites[1])
    mapping = Dict("none"=>-1, "odd"=>1, "even"=>0)
    if parity_symmetry_type == "none" && inversion_symmetry_type == "none"
        return randomMPS(sites; linkdims=Nspec)
    end

    for psi in (randomMPS(sites,[[1 for i in 1:Nsites-1]..., i]) for i=1:Nspec)
        f = flux(psi)
        if parity_symmetry_type != "none" && mapping[parity_symmetry_type] != val(f, "parity")
            continue
        end
        if inversion_symmetry_type != "none" && mapping[inversion_symmetry_type] != val(f, "inv_sym")
            continue
        end
        return psi
    end
end

function parity_symmetry(matrix)
	Nspec, _ = size(matrix)
	is_even = (Nspec+1)%2
	mmax = div(Nspec,2)
	state_to_ind(i) = i + mmax + (1-is_even)
	evens = state_to_ind.(filter(iseven, -mmax + is_even:mmax)) 
	odds = state_to_ind.(filter(isodd, -mmax + is_even:mmax))
	order = [evens..., odds...]
	return matrix[order, order]
end

function m_inversion_symmetry(matrix)
	Nspec, _ = size(matrix)
	@assert Nspec%2 == 1
	operator = zeros(Nspec, Nspec)
	for i in 1:Nspec÷2
		operator[i,i] = 1/sqrt(2)
		operator[i,end-i+1] = 1/sqrt(2)
		operator[end-i+1,end-i+1] = -1/sqrt(2)
		operator[end-i+1,i] = 1/sqrt(2)
	end
	operator[Nspec÷2+1, Nspec÷2+1] = 1
	operator = operator' * matrix *operator
	operator[abs.(operator) .< 1e-10] .= 0
	return operator
end

function dvr_symmetric_basis(matrix)
	# basis which the eigenvalues are rotation and inversion eigenstates
	Nspec, _ = size(matrix)
	mmax = Nspec÷2
	@assert Nspec%2 == 0
	operator = zeros(Nspec, Nspec)
	for (s, ix, iy) in [(1, 0,0),(1, 0,mmax),(1,mmax,0),(-1, mmax,mmax)]
		neg = ((iy == mmax) ? -1 : 1)
		operator[1+ix,1+iy] = s/sqrt(2)
		for i in 1:Nspec÷4
			operator[i+1+ix,i+1+iy] = s/2
			operator[i+1+ix,mmax-i+1+iy] = s/2
			operator[mmax-i+1+ix,mmax-i+1+iy] = -neg*s/2
			operator[mmax-i+1+ix,i+1+iy] = neg*s/2
		end
		if Nspec%4 == 0
			operator[Nspec÷4+1+ix, Nspec÷4+1+iy] = s/sqrt(2)
		end
	end
	operator = operator'*matrix*operator
	operator[abs.(operator) .< 1e-10] .= 0
	return operator
end

function dvr_inversion_symmetry(matrix)
	Nspec, _ = size(matrix)
	operator = zeros(Nspec, Nspec)
	operator[1,1] = 1
	for i in 1:Nspec÷2
		operator[i+1,i+1] = 1/sqrt(2)
		operator[i+1,end-i+1] = 1/sqrt(2)
		operator[end-i+1,end-i+1] = -1/sqrt(2)
		operator[end-i+1,i+1] = 1/sqrt(2)
	end
	if Nspec%2 == 0
		operator[Nspec÷2+1, Nspec÷2+1] = 1
	end
	operator = operator' * matrix *operator
	operator[abs.(operator) .< 1e-10] .= 0
	return operator
end

function dvr_rotation_symmetry(matrix)
	Nspec, _ = size(matrix)
	mmax = Nspec÷2
	@assert Nspec%2 == 0
	operator = zeros(Nspec, Nspec)
	for i in 1:mmax
		operator[i,i] = 1/sqrt(2)
		operator[i,i+mmax] = 1/sqrt(2)
		operator[i+mmax,i] = 1/sqrt(2)
		operator[i+mmax,i+mmax] = -1/sqrt(2)
	end
	operator = operator' * matrix *operator
	operator[abs.(operator) .< 1e-10] .= 0
	return operator
end

function trivial_symmetry(matrix)
	return matrix
end
	