function ITensors.space(::SiteType"PlaRotor";conserve_qns=false) 
    if conserve_qns
        return [QN("T",-1)=>1,QN("T",1)=>1]
    end
    return Nspec
end

#######################################################################
function ITensors.op!(Op::ITensor,::OpName"T",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            Op[s'=>j,s=>i] = T[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"m",::SiteType"PlaRotor" ,s::Index)
	for i=1:Nspec
		for j=1:Nspec
			Op[s'=>j,s=>i] = m[j,i]
		end
	end
end
#######################################################################
function ITensors.op!(Op::ITensor,::OpName"X",::SiteType"PlaRotor" ,s::Index)
	for i=1:Nspec
		for j=1:Nspec
			Op[s'=>j,s=>i] = X[j,i]
		end
	end
end
#######################################################################
function ITensors.op!(Op::ITensor,::OpName"Y",::SiteType"PlaRotor" ,s::Index)
	for i=1:Nspec
		for j=1:Nspec
			Op[s'=>j,s=>i] = Y[j,i]
		end
	end
end
#######################################################################
function ITensors.op!(Op::ITensor,::OpName"Up",::SiteType"PlaRotor" ,s::Index)
	for i=1:Nspec
        for j=1:Nspec
            Op[s'=>j,s=>i] = Up[j,i]
        end
	end     
end
#######################################################################
function ITensors.op!(Op::ITensor,::OpName"Down",::SiteType"PlaRotor" ,s::Index)
	for i=1:Nspec
        for j=1:Nspec
            Op[s'=>j,s=>i] = Down[j,i]
        end
	end             
end
#######################################################################
function ITensors.op!(Op::ITensor,::OpName"Ycomp",::SiteType"PlaRotor" ,s::Index)
	for i=1:Nspec
		for j=1:Nspec
			Op[s'=>j,s=>i] = Y[j,i]*im
		end
	end
end
