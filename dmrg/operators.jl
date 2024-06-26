
##code example for symmetry and variable dim
#=
    function ITensors.space(
  ::SiteType"Qudit";
  dim=2,
  conserve_qns=false,
  conserve_number=conserve_qns,
  qnname_number="Number",
)
  if conserve_number
    return [QN(qnname_number, n - 1) => 1 for n in 1:dim]
  end
  return dim
end
=#

function create_Hamiltonian(g, sites, pairs; Estrength=0,angle=0, evod="m")
	Nsites = length(sites)

	Nsecond = zeros(Int64,(Nsites-1))

	for i=1:Nsites-1
		if pairs == "nearest"
			Nsecond[i]=i+1
		elseif pairs == "allpairs"
			Nsecond[i]=Nsites
        else
            throw(string("Pairs is: ", pairs))
        end
	end

	ampo = AutoMPO()
	for i=1:Nsites-1
		ampo += 1.0,"T",i
		for j=i+1:Nsecond[i]
			c=g/((abs(j-i))^3)
			if evod == "dvr"
				# y_iy_j#
				ampo += 1.0*c,"Y",i,"Y",j
				# 2*x_ix_j#
				ampo += -2.0*c,"X",i,"X",j
			else
				# up up
				ampo +=-.75*c,"Up",i,"Up",j
				# up down
				ampo +=-.25*c,"Up",i,"Down",j
				# down up 
				ampo +=-.25*c,"Down",i,"Up",j
				# down down
				ampo +=-.75*c,"Down",i,"Down",j

                # ampo += -1.0*c,"Y",i,"Y",j
				# ampo += -2.0*c,"X",i,"X",j
			end
		end
		#Electric field#
		
		if !iszero(Estrength)
			ampo += -cos(angle)*Estrength,"X",i
			ampo += -sin(angle)*Estrength*fac2,"Y",i
		end
	end
	ampo += 1.0,"T",Nsites
	#Electric field#

	if !iszero(Estrength)
		ampo += -cos(angle)*Estrength,"X",Nsites
		ampo += -sin(angle)*Estrength*fac2,"Y",Nsites
	end
	H = MPO(ampo,sites)
	return H
end

function label_states_by_parity(dim::Int)
    is_even = (dim+1)%2
    even = []
    odd = []
    k=0
    mmax=div(dim,2)
    for m=-mmax+is_even:mmax
        k+=1
        if mod(abs(m),2) == 0
            append!(even,k)
        else
            append!(odd,k)
        end
    end
    return odd,even
end

function ITensors.space(
    ::SiteType"PlaRotor";
    dim,
    basis="m",
    conserve_parity=false,
    conserve_inversion_symmetry=false,
    )
    mmax = dim ÷ 2
    is_odd = dim%2
    is_even = 1 - is_odd

    # no conserved quantities
    if !conserve_parity && !conserve_inversion_symmetry
        return dim
    end

    if basis == "m"
        if conserve_inversion_symmetry && conserve_parity
            @assert is_odd == 1
            quantum_numbers = Pair{QN, Int64}[]
            for m = -mmax:mmax
                push!(quantum_numbers, QN(("parity", m%2,2), ("inv_sym", Int(m>0), 2))=>1)
            end
            return quantum_numbers
        elseif conserve_inversion_symmetry
            @assert is_odd == 1
            return [QN("inv_sym",0,2)=>mmax+1, QN("inv_sym",1,2)=>mmax]
        else
            quantum_numbers = Pair{QN, Int64}[]
            for m = -mmax+is_even:mmax
                push!(quantum_numbers, QN("parity", m%2, 2)=>1)
            end
            return quantum_numbers
        end
    elseif basis == "dvr"
        if conserve_inversion_symmetry && conserve_parity
            @assert is_even == 1
            return [
                QN(("parity",0,2), ("inv_sym",0,2))=>1 + mmax÷2, 
                QN(("parity",0,2), ("inv_sym",1,2))=>(mmax-1)÷2, 
                QN(("parity",1,2), ("inv_sym",0,2))=>(mmax+1)÷2, 
                QN(("parity",1,2), ("inv_sym",1,2))=>mmax÷2, 
            ]
        elseif conserve_inversion_symmetry
            return [QN("inv_sym",0,2)=>mmax+1, QN("inv_sym",1,2)=>mmax-is_even]
        else
            @assert is_even == 1
            return [QN("parity",0,2)=>mmax, QN("parity",1,2)=>mmax]
        end
    else
        throw("Invalid basis: $basis. Must be m or dvr")
    end
end

#######################################################################
function ITensors.op!(Op::ITensor,::OpName"T",::SiteType"PlaRotor" ,s::Index)
#@show T
    for i=1:Nspec
        for j=1:Nspec
            iszero(T[j,i]) ? nothing : Op[s'=>j,s=>i] = T[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"m",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(m[j,i]) ? nothing : Op[s'=>j,s=>i] = m[j,i]
        end
    end
end
#######################################################################
function ITensors.op!(Op::ITensor,::OpName"X",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(X[j,i]) ? nothing : Op[s'=>j,s=>i] = X[j,i]
        end
    end
end
#######################################################################
function ITensors.op!(Op::ITensor,::OpName"Y",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(Y[j,i]) ? nothing : Op[s'=>j,s=>i] = Y[j,i]
        end
    end
end
#######################################################################
function ITensors.op!(Op::ITensor,::OpName"Up",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(Up[j,i]) ? nothing : Op[s'=>j,s=>i] = Up[j,i]
        end
    end     
end
#######################################################################
function ITensors.op!(Op::ITensor,::OpName"Down",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(Down[j,i]) ? nothing : Op[s'=>j,s=>i] = Down[j,i]
        end
    end         
end
# #######################################################################
# function ITensors.op!(Op::ITensor,::OpName"Ycomp",::SiteType"PlaRotor" ,s::Index)
#     for i=1:Nspec
#         for j=1:Nspec
#             iszero(Y[j,i]) ? nothing : Op[s'=>j,s=>i] = Y[j,i]*im
#         end
#     end
# end

function ITensors.op!(Op::ITensor,::OpName"mInvert",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(mInvert[j,i]) ? nothing : Op[s'=>j,s=>i] = mInvert[j,i]
        end
    end
end


function ITensors.op!(Op::ITensor,::OpName"SmallEProj",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(SmallEProj[j,i]) ? nothing : Op[s'=>j,s=>i] = SmallEProj[j,i]
        end
    end
end

