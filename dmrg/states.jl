#module states
#export symmetry
################################################
#=
function symmetry(matrix,state,mmax)
even = []
odd = []
k=0
for m=-mmax:mmax
	k+=1
	if mod(abs(m),2) == 0
		append!(even,k)
	else
		append!(odd,k)
	end
end
indices = []
if state == "even"
	indices = sort(even)
elseif state == "odd"
	indices = sort(odd)
end
Nspec = length(indices)
matrix_new = zeros(Nspec,Nspec)
i=0
for i1 in indices
	i+=1
	j=0
	for j1 in indices
			j+=1
			matrix_new[i,j]=matrix[i1,j1]
	end
end
return matrix_new
end
=#
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

function trivial_symmetry(matrix)
	return matrix
end
	


################################################
#end
