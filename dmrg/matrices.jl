#module matrices

using LinearAlgebra

#export kinetic,Xoperator,Yoperator,Upoperator,Downoperator,PNoperator,m
#################################################################
function kinetic(mmax)

	matrix = zeros((2*mmax+1),(2*mmax+1))

	k=0
	for m=-mmax:mmax
		k+=1
		matrix[k,k] = m*m
	end	

	return matrix
end
function m(mmax)
	matrix = zeros((2*mmax+1),(2*mmax+1))
	k=0
	for m=-mmax:mmax
		k+=1
		matrix[k,k] = m
	end	
	return matrix
end
#################################################################
function Xoperator(mmax)

	diagonal = zeros((2*mmax+1))
	subdiagonal = 0.5*ones((2*mmax))

	matrix=Tridiagonal(subdiagonal,diagonal,subdiagonal)

	return matrix
end
#################################################################
function Yoperator(mmax)

	diagonal = zeros((2*mmax+1))
	upper_subdiagonal = 0.5*ones((2*mmax))
	lower_subdiagonal = -0.5*ones((2*mmax))

	matrix=Tridiagonal(lower_subdiagonal,diagonal,upper_subdiagonal)

	return matrix
end
#################################################################
function Upoperator(mmax)

	diagonal = zeros((2*mmax+1))
	upper_subdiagonal = zeros((2*mmax))
	subdiagonal = ones((2*mmax))
		
	matrix=Tridiagonal(subdiagonal,diagonal,upper_subdiagonal)
		
	return matrix
end
#################################################################
function Downoperator(mmax)

	diagonal = zeros((2*mmax+1))
	upper_subdiagonal = ones((2*mmax))
	subdiagonal = zeros((2*mmax))

	matrix=Tridiagonal(subdiagonal,diagonal,upper_subdiagonal)
		
	return matrix
end
function MInversionOperator(mmax)
	return rotl90(Diagonal(ones(2*mmax+1)))
end
# #################################################################
# function InversionProjector(mmax; symmetry="even")
# 	diagonal = Diagonal(ones(2*mmax+1)./2)
# 	if symmetry == "even"
# 		return diagonal + rotl90(diagonal)
# 	else
# 		return diagonal - rotl90(diagonal)
# 	end
# end

#################################################################
function SmallEnergyProjector(mmax; m=1)
	diagonal = zeros(2*mmax+1)
	diagonal[(mmax-m+1):(mmax+m+1)] = ones(2*m+1)
	return Diagonal(diagonal)
end

function ReflectionOperator(mmax)
	dim = (2*mmax+1)
	arr = zeros(dim,dim,dim,dim)
	for i = 1:dim
		for j = 1:dim
			arr[i,j,j,i] = 1
		end
	end
	return reshape(arr, dim*dim,dim*dim)
end

#end 
