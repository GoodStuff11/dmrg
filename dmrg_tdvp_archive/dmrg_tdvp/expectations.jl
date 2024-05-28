function vN_entropy(wf,mbond)
	if length(wf) == 2
		orthogonalize!(wf, 2)
		U,S,V = svd(wf[2], (siteind(wf,2)))
	else
		orthogonalize!(wf, mbond)
		U,S,V = svd(wf[mbond], (linkind(wf, mbond-1), siteind(wf,mbond)))
	end
	SvN = 0.0
	renyi = 0.0
	schmidtvalues=zeros(dim(S, 1))
	for n=1:dim(S, 1)
		p = S[n,n]^2
		schmidtvalues[n]=p
		SvN -= p * log(p)
		renyi += p^2
	end
	renyi=-0.5*log(renyi)

	return SvN,renyi,schmidtvalues
end
################################################################################
function polarization(wf,Nsites,evod)
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
function correlation(wf,Nsites,evod)
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