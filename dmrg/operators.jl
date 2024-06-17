
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


function create_Hamiltonian(g, sites, Nsecond; Estrength=0,angle=0, fac1=1, fac2=1)
	Nsites = length(sites)
	ampo = AutoMPO()
	for i=1:Nsites-1
		ampo += 1.0,"T",i
		for j=i+1:Nsecond[i]
			c=g/((abs(j-i))^3)
			if evod == "dvr"
				# y_iy_j#
				ampo += 1.0*c*fac1,"Y",i,"Y",j
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

function ITensors.space(
    ::SiteType"PlaRotor";
    dim=3,
    conserve_qns=false,
    conserve_parity=false,
    conserve_L=false,
    qnname_parity="Parity",
    qnname_totalL ="L"
    )
    ##currently only parity supported
    ##
    if conserve_qns
        conserve_parity=true
        conserve_l=true
    end
    if conserve_parity || conserve_L
        mmax=div(dim,2)
        
        #evenstates,oddstates=symmetry(dim)
        #this requires reordering of states in the definitions
        [QN(qnname_parity,0,2)=>length(filter(iseven,-mmax:mmax)),QN(qnname_parity,1,2)=>length(filter(isodd,-mmax:mmax))]
        #this does not but leads to fragmented blocks
        #return [QN(qnname_parity,Int(isodd(i)),2)=>1 for i in -mmax:mmax]
    elseif conserve_parity && conserve_L
        mmax=div(dim,2)
        [QN((qnname_parity,isodd(m),2),(qnname_totalL,m,1))=>1 for m in -mmax:mmax]
    elseif conserve_L
        [QN(qnname_totalL,m,1)=>1 for m in -mmax:mmax]
    else
        return dim
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



### GENERATED CODE ###
function ITensors.op!(Op::ITensor,::OpName"U_0_0",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_0_0[j,i]) ? nothing : Op[s'=>j,s=>i] = U_0_0[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_0_1",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_0_1[j,i]) ? nothing : Op[s'=>j,s=>i] = U_0_1[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_0_2",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_0_2[j,i]) ? nothing : Op[s'=>j,s=>i] = U_0_2[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_0_3",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_0_3[j,i]) ? nothing : Op[s'=>j,s=>i] = U_0_3[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_0_4",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_0_4[j,i]) ? nothing : Op[s'=>j,s=>i] = U_0_4[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_0_5",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_0_5[j,i]) ? nothing : Op[s'=>j,s=>i] = U_0_5[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_0_6",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_0_6[j,i]) ? nothing : Op[s'=>j,s=>i] = U_0_6[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_0_7",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_0_7[j,i]) ? nothing : Op[s'=>j,s=>i] = U_0_7[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_0_8",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_0_8[j,i]) ? nothing : Op[s'=>j,s=>i] = U_0_8[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_0_9",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_0_9[j,i]) ? nothing : Op[s'=>j,s=>i] = U_0_9[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_0_10",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_0_10[j,i]) ? nothing : Op[s'=>j,s=>i] = U_0_10[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_1_0",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_1_0[j,i]) ? nothing : Op[s'=>j,s=>i] = U_1_0[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_1_1",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_1_1[j,i]) ? nothing : Op[s'=>j,s=>i] = U_1_1[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_1_2",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_1_2[j,i]) ? nothing : Op[s'=>j,s=>i] = U_1_2[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_1_3",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_1_3[j,i]) ? nothing : Op[s'=>j,s=>i] = U_1_3[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_1_4",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_1_4[j,i]) ? nothing : Op[s'=>j,s=>i] = U_1_4[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_1_5",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_1_5[j,i]) ? nothing : Op[s'=>j,s=>i] = U_1_5[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_1_6",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_1_6[j,i]) ? nothing : Op[s'=>j,s=>i] = U_1_6[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_1_7",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_1_7[j,i]) ? nothing : Op[s'=>j,s=>i] = U_1_7[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_1_8",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_1_8[j,i]) ? nothing : Op[s'=>j,s=>i] = U_1_8[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_1_9",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_1_9[j,i]) ? nothing : Op[s'=>j,s=>i] = U_1_9[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_1_10",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_1_10[j,i]) ? nothing : Op[s'=>j,s=>i] = U_1_10[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_2_0",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_2_0[j,i]) ? nothing : Op[s'=>j,s=>i] = U_2_0[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_2_1",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_2_1[j,i]) ? nothing : Op[s'=>j,s=>i] = U_2_1[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_2_2",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_2_2[j,i]) ? nothing : Op[s'=>j,s=>i] = U_2_2[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_2_3",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_2_3[j,i]) ? nothing : Op[s'=>j,s=>i] = U_2_3[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_2_4",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_2_4[j,i]) ? nothing : Op[s'=>j,s=>i] = U_2_4[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_2_5",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_2_5[j,i]) ? nothing : Op[s'=>j,s=>i] = U_2_5[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_2_6",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_2_6[j,i]) ? nothing : Op[s'=>j,s=>i] = U_2_6[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_2_7",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_2_7[j,i]) ? nothing : Op[s'=>j,s=>i] = U_2_7[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_2_8",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_2_8[j,i]) ? nothing : Op[s'=>j,s=>i] = U_2_8[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_2_9",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_2_9[j,i]) ? nothing : Op[s'=>j,s=>i] = U_2_9[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_2_10",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_2_10[j,i]) ? nothing : Op[s'=>j,s=>i] = U_2_10[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_3_0",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_3_0[j,i]) ? nothing : Op[s'=>j,s=>i] = U_3_0[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_3_1",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_3_1[j,i]) ? nothing : Op[s'=>j,s=>i] = U_3_1[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_3_2",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_3_2[j,i]) ? nothing : Op[s'=>j,s=>i] = U_3_2[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_3_3",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_3_3[j,i]) ? nothing : Op[s'=>j,s=>i] = U_3_3[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_3_4",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_3_4[j,i]) ? nothing : Op[s'=>j,s=>i] = U_3_4[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_3_5",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_3_5[j,i]) ? nothing : Op[s'=>j,s=>i] = U_3_5[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_3_6",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_3_6[j,i]) ? nothing : Op[s'=>j,s=>i] = U_3_6[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_3_7",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_3_7[j,i]) ? nothing : Op[s'=>j,s=>i] = U_3_7[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_3_8",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_3_8[j,i]) ? nothing : Op[s'=>j,s=>i] = U_3_8[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_3_9",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_3_9[j,i]) ? nothing : Op[s'=>j,s=>i] = U_3_9[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_3_10",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_3_10[j,i]) ? nothing : Op[s'=>j,s=>i] = U_3_10[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_4_0",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_4_0[j,i]) ? nothing : Op[s'=>j,s=>i] = U_4_0[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_4_1",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_4_1[j,i]) ? nothing : Op[s'=>j,s=>i] = U_4_1[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_4_2",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_4_2[j,i]) ? nothing : Op[s'=>j,s=>i] = U_4_2[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_4_3",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_4_3[j,i]) ? nothing : Op[s'=>j,s=>i] = U_4_3[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_4_4",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_4_4[j,i]) ? nothing : Op[s'=>j,s=>i] = U_4_4[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_4_5",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_4_5[j,i]) ? nothing : Op[s'=>j,s=>i] = U_4_5[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_4_6",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_4_6[j,i]) ? nothing : Op[s'=>j,s=>i] = U_4_6[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_4_7",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_4_7[j,i]) ? nothing : Op[s'=>j,s=>i] = U_4_7[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_4_8",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_4_8[j,i]) ? nothing : Op[s'=>j,s=>i] = U_4_8[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_4_9",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_4_9[j,i]) ? nothing : Op[s'=>j,s=>i] = U_4_9[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_4_10",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_4_10[j,i]) ? nothing : Op[s'=>j,s=>i] = U_4_10[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_5_0",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_5_0[j,i]) ? nothing : Op[s'=>j,s=>i] = U_5_0[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_5_1",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_5_1[j,i]) ? nothing : Op[s'=>j,s=>i] = U_5_1[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_5_2",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_5_2[j,i]) ? nothing : Op[s'=>j,s=>i] = U_5_2[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_5_3",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_5_3[j,i]) ? nothing : Op[s'=>j,s=>i] = U_5_3[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_5_4",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_5_4[j,i]) ? nothing : Op[s'=>j,s=>i] = U_5_4[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_5_5",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_5_5[j,i]) ? nothing : Op[s'=>j,s=>i] = U_5_5[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_5_6",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_5_6[j,i]) ? nothing : Op[s'=>j,s=>i] = U_5_6[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_5_7",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_5_7[j,i]) ? nothing : Op[s'=>j,s=>i] = U_5_7[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_5_8",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_5_8[j,i]) ? nothing : Op[s'=>j,s=>i] = U_5_8[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_5_9",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_5_9[j,i]) ? nothing : Op[s'=>j,s=>i] = U_5_9[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_5_10",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_5_10[j,i]) ? nothing : Op[s'=>j,s=>i] = U_5_10[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_6_0",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_6_0[j,i]) ? nothing : Op[s'=>j,s=>i] = U_6_0[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_6_1",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_6_1[j,i]) ? nothing : Op[s'=>j,s=>i] = U_6_1[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_6_2",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_6_2[j,i]) ? nothing : Op[s'=>j,s=>i] = U_6_2[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_6_3",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_6_3[j,i]) ? nothing : Op[s'=>j,s=>i] = U_6_3[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_6_4",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_6_4[j,i]) ? nothing : Op[s'=>j,s=>i] = U_6_4[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_6_5",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_6_5[j,i]) ? nothing : Op[s'=>j,s=>i] = U_6_5[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_6_6",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_6_6[j,i]) ? nothing : Op[s'=>j,s=>i] = U_6_6[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_6_7",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_6_7[j,i]) ? nothing : Op[s'=>j,s=>i] = U_6_7[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_6_8",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_6_8[j,i]) ? nothing : Op[s'=>j,s=>i] = U_6_8[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_6_9",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_6_9[j,i]) ? nothing : Op[s'=>j,s=>i] = U_6_9[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_6_10",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_6_10[j,i]) ? nothing : Op[s'=>j,s=>i] = U_6_10[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_7_0",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_7_0[j,i]) ? nothing : Op[s'=>j,s=>i] = U_7_0[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_7_1",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_7_1[j,i]) ? nothing : Op[s'=>j,s=>i] = U_7_1[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_7_2",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_7_2[j,i]) ? nothing : Op[s'=>j,s=>i] = U_7_2[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_7_3",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_7_3[j,i]) ? nothing : Op[s'=>j,s=>i] = U_7_3[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_7_4",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_7_4[j,i]) ? nothing : Op[s'=>j,s=>i] = U_7_4[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_7_5",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_7_5[j,i]) ? nothing : Op[s'=>j,s=>i] = U_7_5[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_7_6",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_7_6[j,i]) ? nothing : Op[s'=>j,s=>i] = U_7_6[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_7_7",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_7_7[j,i]) ? nothing : Op[s'=>j,s=>i] = U_7_7[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_7_8",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_7_8[j,i]) ? nothing : Op[s'=>j,s=>i] = U_7_8[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_7_9",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_7_9[j,i]) ? nothing : Op[s'=>j,s=>i] = U_7_9[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_7_10",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_7_10[j,i]) ? nothing : Op[s'=>j,s=>i] = U_7_10[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_8_0",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_8_0[j,i]) ? nothing : Op[s'=>j,s=>i] = U_8_0[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_8_1",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_8_1[j,i]) ? nothing : Op[s'=>j,s=>i] = U_8_1[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_8_2",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_8_2[j,i]) ? nothing : Op[s'=>j,s=>i] = U_8_2[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_8_3",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_8_3[j,i]) ? nothing : Op[s'=>j,s=>i] = U_8_3[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_8_4",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_8_4[j,i]) ? nothing : Op[s'=>j,s=>i] = U_8_4[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_8_5",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_8_5[j,i]) ? nothing : Op[s'=>j,s=>i] = U_8_5[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_8_6",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_8_6[j,i]) ? nothing : Op[s'=>j,s=>i] = U_8_6[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_8_7",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_8_7[j,i]) ? nothing : Op[s'=>j,s=>i] = U_8_7[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_8_8",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_8_8[j,i]) ? nothing : Op[s'=>j,s=>i] = U_8_8[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_8_9",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_8_9[j,i]) ? nothing : Op[s'=>j,s=>i] = U_8_9[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_8_10",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_8_10[j,i]) ? nothing : Op[s'=>j,s=>i] = U_8_10[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_9_0",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_9_0[j,i]) ? nothing : Op[s'=>j,s=>i] = U_9_0[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_9_1",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_9_1[j,i]) ? nothing : Op[s'=>j,s=>i] = U_9_1[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_9_2",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_9_2[j,i]) ? nothing : Op[s'=>j,s=>i] = U_9_2[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_9_3",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_9_3[j,i]) ? nothing : Op[s'=>j,s=>i] = U_9_3[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_9_4",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_9_4[j,i]) ? nothing : Op[s'=>j,s=>i] = U_9_4[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_9_5",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_9_5[j,i]) ? nothing : Op[s'=>j,s=>i] = U_9_5[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_9_6",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_9_6[j,i]) ? nothing : Op[s'=>j,s=>i] = U_9_6[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_9_7",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_9_7[j,i]) ? nothing : Op[s'=>j,s=>i] = U_9_7[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_9_8",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_9_8[j,i]) ? nothing : Op[s'=>j,s=>i] = U_9_8[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_9_9",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_9_9[j,i]) ? nothing : Op[s'=>j,s=>i] = U_9_9[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_9_10",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_9_10[j,i]) ? nothing : Op[s'=>j,s=>i] = U_9_10[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_10_0",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_10_0[j,i]) ? nothing : Op[s'=>j,s=>i] = U_10_0[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_10_1",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_10_1[j,i]) ? nothing : Op[s'=>j,s=>i] = U_10_1[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_10_2",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_10_2[j,i]) ? nothing : Op[s'=>j,s=>i] = U_10_2[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_10_3",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_10_3[j,i]) ? nothing : Op[s'=>j,s=>i] = U_10_3[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_10_4",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_10_4[j,i]) ? nothing : Op[s'=>j,s=>i] = U_10_4[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_10_5",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_10_5[j,i]) ? nothing : Op[s'=>j,s=>i] = U_10_5[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_10_6",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_10_6[j,i]) ? nothing : Op[s'=>j,s=>i] = U_10_6[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_10_7",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_10_7[j,i]) ? nothing : Op[s'=>j,s=>i] = U_10_7[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_10_8",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_10_8[j,i]) ? nothing : Op[s'=>j,s=>i] = U_10_8[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_10_9",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_10_9[j,i]) ? nothing : Op[s'=>j,s=>i] = U_10_9[j,i]
        end
    end
end
function ITensors.op!(Op::ITensor,::OpName"U_10_10",::SiteType"PlaRotor" ,s::Index)
    for i=1:Nspec
        for j=1:Nspec
            iszero(U_10_10[j,i]) ? nothing : Op[s'=>j,s=>i] = U_10_10[j,i]
        end
    end
end
### END OF GENERATED CODE ###

