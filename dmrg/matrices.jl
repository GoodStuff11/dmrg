#module matrices

using LinearAlgebra

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
	@assert dim % 2 == 1
	# this symmetry doesn't work for even dim
	# works 
	return rotl90(Diagonal(ones(dim)))
end

function MParityOperator(dim)
	@assert dim % 2 == 1
	mmax = dim ÷ 2
	arr = zeros(dim,dim)
	for i = 1:dim
		arr[i, i] = (-1)^((i-mmax)%2)
	end
	return arr
end

function phiReflectionOperator(dim)
	# phi -> -phi
	arr = zeros(dim,dim)
	arr[1,1] = 1
	arr[2:end,2:end] = rotl90(Diagonal(ones(dim-1)))
	return arr
end


function phiRotationOperator(dim)
	@assert dim % 2 == 0
	arr = zeros(dim,dim)
	for i = 1:dim
		arr[i, (i+dim÷2-1)%dim+1] = 1
	end
	return arr
end

function SmallEnergyProjector(dim; m=1)
	# only for m basis
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
	return reshape(arr,dim*dim,dim*dim)
end

#end 
