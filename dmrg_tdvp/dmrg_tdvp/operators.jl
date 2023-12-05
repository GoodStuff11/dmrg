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

function label_states_by_parity(dim::Int)
    even = []
    odd = []
    k=0
    mmax=div(dim,2)
    for m=-mmax:mmax
        k+=1
        if mod(abs(m),2) == 0
            append!(even,k)
        else
            append!(odd,k)
        end
    end
    return odd,even
end

function parity_symmetry(matrix)
    Nspec,_ = size(matrix)
    mmax = div(Nspec,2)
    state_to_ind(i)=i+mmax+1
    evens=state_to_ind.(filter(iseven,-mmax:mmax)) 
    odds=state_to_ind.(filter(isodd,-mmax:mmax))
    order=[evens...,odds...]
    matrix_new = zeros(Nspec,Nspec)
    return matrix[order,order]
end

function trivial_symmetry(matrix)
    return matrix
end

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

function set_operators(dim::Int64; evod::String="all")
    global Nspec = dim
    mmax = div(dim,2)
    
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

function ITensors.space(::SiteType"PlaRotor";
    dim=3,
    conserve_qns=false,
    conserve_parity=false,
    conserve_L=false,
    qnname_parity="Parity",
    qnname_totalL ="L"
    )

    if !isnothing(dim)
        set_operators(dim)
    elseif isnothing(Nspec)
        error("Either input dim=mmax indo siteinds or call set_operators(mmax)")
    end
    if conserve_qns
        conserve_parity = true
        conserve_L = true
    end
    if conserve_parity && conserve_L
        mmax=div(dim,2)
        return [QN((qnname_parity,isodd(m),2),(qnname_totalL,m,1))=>1 for m in -mmax:mmax]
    elseif conserve_parity
        mmax=div(dim,2)
        
        #evenstates,oddstates=symmetry(dim)
        #this requires reordering of states in the definitions
        return [QN(qnname_parity,0,2)=>length(filter(iseven,-mmax:mmax)),QN(qnname_parity,1,2)=>length(filter(isodd,-mmax:mmax))]
        #this does not but leads to fragmented blocks
        #return [QN(qnname_parity,Int(isodd(i)),2)=>1 for i in -mmax:mmax]
    elseif conserve_L
        return [QN(qnname_totalL,m,1)=>1 for m in -mmax:mmax]
    else
        return dim
    end
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