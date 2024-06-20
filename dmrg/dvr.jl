#module dvr

#export exp_dvr
#############################################################################

function exp_dvr(dim)
	N = dim÷2
	#Write DVR grid and position operators on grid#
	phi = zeros(dim)
	X = zeros((dim,dim))
	Y = zeros((dim,dim))
	
	for ii=1:dim
		phi[ii] = ii*2.0*pi/dim
		X[ii,ii] = cos(phi[ii])	
		Y[ii,ii] = sin(phi[ii])	
	end
	
	#Write kinetic energy matrix in DVR basis#
	T = zeros((dim,dim))
	
	for ii=1:dim
		for jj=ii:dim
			if dim%2 == 0
				T[ii,jj] = 1/(2.0*sin(pi*(ii-jj)/dim)^2)*(-1.)^(ii-jj)
			else
				T[ii,jj] = cos(pi*(ii-jj)/dim)/(2.0*sin(pi*(ii-jj)/dim)^2)*(-1.)^(ii-jj)
			end
			T[jj,ii] = T[ii,jj]
		end
		if dim%2 == 0
			T[ii,ii] = (2N^2+1)/6
		else
			T[ii,ii] = N*(N+1)/3
		end
	end
	
	return T,X,Y
end
#############################################################################
#end
