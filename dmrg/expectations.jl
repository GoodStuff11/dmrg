#module expectations
  
#using ITensors
#using LinearAlgebra


#export vN_entropy,polarization,correlation
################################################################################
function vN_entropy(_wf,mbond)
	wf=copy(_wf)	###to ensure that this function is not mutating wf
	if length(wf) == 2
		orthogonalize!(wf, 2)
		U,S,V = svd(wf[2], (siteind(wf,2)))
	else
		orthogonalize!(wf, mbond)
		U,S,V = svd(wf[mbond], (linkind(wf, mbond-1), siteind(wf,mbond)))
	end
	SvN = 0.0
	purity = 0.0
	schmidtvalues=zeros(dim(S, 1))
	for n=1:dim(S, 1)
		schmidtvalues[n] = S[n,n]^2
	end
	# divide by the sum of probabilities to ensure probabilities are normalized
	sum_probability = sum(schmidtvalues)
	for n=1:dim(S, 1)
		schmidtvalues[n] /= sum_probability
		p = schmidtvalues[n]
		SvN -= p>0 ? p * log(p) : 0
		purity += p^2
	end
	renyi=-0.5*log(purity)
	return SvN,purity, schmidtvalues
	end
################################################################################
function polarization(wf,Nsites,Nbasis,evod,Xmat,Ymat)

	global X = Xmat
	global Y = Ymat
	global Nspec = Nbasis

	include("operators.jl")

	if evod == "dvr"
		mux = expect(wf,"X")
		muy = expect(wf,"Y")

		dumX = correlation_matrix(wf,"X","X")
		dumY = correlation_matrix(wf,"Y","Y")
	else
		mux = expect(wf,"X")
		#muy = expect(wf,"Ycomp")
		muy = expect(wf,"Y")
		
		dumX = correlation_matrix(wf,"X","X")
		#dumY = correlation_matrix(wf,"Ycomp","Ycomp")
		dumY = correlation_matrix(wf,"Y","Y")
	end

	mux2=0.0
	muy2=0.0
	for ii=1:Nsites
		mux2+=dumX[ii,ii]
		if evod == "dvr"
			muy2+=dumY[ii,ii]
		else
			muy2-=dumY[ii,ii]
		end
	end

return real(sum(mux)),real(sum(muy))
end
################################################################################
function correlation(wf,Nsites,Nbasis,evod,Xmat,Ymat)

	global X = Xmat
	global Y = Ymat
	global Nspec = Nbasis

	include("operators.jl")

	if evod == "dvr"
		dumX = correlation_matrix(wf,"X","X")
		dumY = correlation_matrix(wf,"Y","Y")
	else
		dumX = correlation_matrix(wf,"X","X")
		#dumY = correlation_matrix(wf,"Ycomp","Ycomp")
		#pn hack
		dumY = correlation_matrix(wf,"Y","Y")
	end

	Xcorr=0.0
	Ycorr=0.0
	for ii=1:Nsites-1
		Xcorr+=dumX[ii,ii+1]
		if evod == "dvr"
			Ycorr+=dumY[ii,ii+1]
		else
			Ycorr-=dumY[ii,ii+1]
		end

	end

return real(Xcorr),real(Ycorr)
end
################################################################################
#end
