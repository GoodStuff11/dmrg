module states
export symmetry
################################################
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
################################################
end
