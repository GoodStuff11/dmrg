module matrices

using LinearAlgebra

export kinetic,Xoperator,Yoperator,Upoperator,Downoperator,PNoperator,m
#################################################################
"""
K|-3>=9|-3>

Creates matrix of form

9000000
0400000
0010000
0000000
0000100
0000040
0000009
"""
function kinetic(mmax)

	matrix = zeros((2*mmax+1),(2*mmax+1))

	k=0
	for m=-mmax:mmax
		k+=1
		matrix[k,k] = m*m
	end	

	return matrix
end

"""
m|-3>=-3|-3>

Creates matrix of form

-3000000
0-200000
00-10000
0000000
0000100
0000040
0000009
"""
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
"""


Creates matrix of form

0 0.5 0 0 0
0.5 0 0.5 0 0
0 0.5 0 0.5 0
0 0 0.5 0 0.5
0 0 0 0.5 0

"""
function Xoperator(mmax)

	diagonal = zeros((2*mmax+1))
	subdiagonal = 0.5*ones((2*mmax))

	matrix=Tridiagonal(subdiagonal,diagonal,subdiagonal)

	return matrix
end
#################################################################
"""

Creates matrix of form

0 0.5 0 0 0
-0.5 0 0.5 0 0
0 -0.5 0 0.5 0
0 0 -0.5 0 0.5
0 0 0 -0.5 0

"""
function Yoperator(mmax)

	diagonal = zeros((2*mmax+1))
	upper_subdiagonal = 0.5*ones((2*mmax))
	lower_subdiagonal = -0.5*ones((2*mmax))

	matrix=Tridiagonal(lower_subdiagonal,diagonal,upper_subdiagonal)

	return matrix
end
#################################################################
"""
U|i>=|i+1>

Creates matrix of form

0 0 0 0 0
1 0 0 0 0
0 1 0 0 0
0 0 1 0 0
0 0 0 1 0

"""
function Upoperator(mmax)

	diagonal = zeros((2*mmax+1))
	upper_subdiagonal = zeros((2*mmax))
	subdiagonal = ones((2*mmax))
		
	matrix=Tridiagonal(subdiagonal,diagonal,upper_subdiagonal)
		
	return matrix
	end
#################################################################
"""
U|i>=|i-1>

Creates matrix of form

0 1 0 0 0
0 0 1 0 0
0 0 0 1 0
0 0 0 0 1
0 0 0 0 0

"""
function Downoperator(mmax)

	diagonal = zeros((2*mmax+1))
	upper_subdiagonal = ones((2*mmax))
	subdiagonal = zeros((2*mmax))

	matrix=Tridiagonal(subdiagonal,diagonal,upper_subdiagonal)
		
	return matrix
end

#################################################################




end 
