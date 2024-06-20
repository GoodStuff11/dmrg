#module matrices

using LinearAlgebra
using SparseArrays

#export kinetic,Xoperator,Yoperator,Upoperator,Downoperator,PNoperator,m
#################################################################
function kinetic(dim)
	matrix = zeros(dim, dim)

	k=0
	for m=-(dim-1)÷2:dim÷2
		k+=1
		matrix[k,k] = m*m
	end	

	return matrix
end
function m(dim)
	matrix = zeros(dim, dim)
	k=0
	for m=-(dim-1)÷2:dim÷2
		k+=1
		matrix[k,k] = m
	end	
	return matrix
end
#################################################################
function Xoperator(dim)

	diagonal = zeros(dim)
	subdiagonal = 0.5*ones(dim - 1)

	matrix=Tridiagonal(subdiagonal,diagonal,subdiagonal)

	return matrix
end
#################################################################
function Yoperator(dim)

	diagonal = zeros(dim)
	upper_subdiagonal = 0.5*ones(dim-1)

	matrix=Tridiagonal(-upper_subdiagonal,diagonal,upper_subdiagonal)

	return matrix
end
#################################################################
function Upoperator(dim)

	diagonal = zeros(dim)
	upper_subdiagonal = zeros(dim-1)
	subdiagonal = ones(dim-1)
		
	matrix=Tridiagonal(subdiagonal,diagonal,upper_subdiagonal)
		
	return matrix
end
#################################################################
function Downoperator(dim)

	diagonal = zeros(dim)
	upper_subdiagonal = ones(dim-1)
	subdiagonal = zeros(dim-1)

	matrix=Tridiagonal(subdiagonal,diagonal,upper_subdiagonal)
		
	return matrix
end
function MInversionOperator(dim)
	# this symmetry doesn't work for even dim
	return rotl90(Diagonal(ones(dim)))
end

function SmallEnergyProjector(dim; m=1)
	mmax = dim÷2
	diagonal = zeros(dim)
	diagonal[(mmax-m+(dim%2)):(mmax+m+(dim%2))] .= 1
	return Diagonal(diagonal)
end

function ReflectionOperator(dim)
	arr = zeros(dim,dim,dim,dim)
	for i = 1:dim
		for j = 1:dim
			arr[i,j,j,i] = 1
		end
	end
	return sparse(reshape(arr, dim*dim,dim*dim))
end

#end 
