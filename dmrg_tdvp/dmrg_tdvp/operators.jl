module operators

using ITensors
using Printf

include("observer.jl")
include("matrices.jl")
include("dvr.jl")
include("expectations.jl")

using .matrices
using .dvr

export space, op!, checkdone!, DemoObserver, vN_entropy, polarization, correlation

Nspec = nothing
T = nothing
m = nothing
X = nothing
Y = nothing
Up = nothing
Down = nothing

function set_operators(mmax::Int64; evod::String="all")
	global Nspec = 2*mmax+1
	Ttmp = kinetic(mmax)
    Xtmp = Xoperator(mmax)
    Ytmp = Yoperator(mmax)
    Uptmp = Upoperator(mmax)
    Downtmp = Downoperator(mmax)

	if evod == "all"
        global T = Ttmp
        global X = Xtmp
        global Y = Ytmp
        global Up = Uptmp
        global Down = Downtmp
    elseif evod == "dvr"
        tmp1,tmp2,tmp3 = exp_dvr(mmax)
        global T = tmp1
        global X = tmp2
        global Y = tmp3
    else
        global T = symmetry(Ttmp,evod,mmax)
        global X = symmetry(Xtmp,evod,mmax)
        global Y = symmetry(Ytmp,evod,mmax)
    end
end

function ITensors.space(::SiteType"PlaRotor";conserve_qns=false, dim=nothing)
	if !isnothing(dim)
		set_operators(dim)
	elseif isnothing(Nspec)
		error("Either input dim=mmax indo siteinds or call set_operators(mmax)")
	end
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
end