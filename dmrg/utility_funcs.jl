#module utility_funcs

using Printf

#export create_file, write_output, write_text

function create_file(path,evod,mmax,Nsites)
    f=open(path,"w")
    if evod == "dvr"
        println(f,"# Ngrid= ",mmax)
        println(f,"# DVR basis")
    else
        println(f,"# mmax= ",mmax)
        println(f,"# m-states: ",evod," states (FBR)")
    end
    println(f,"# Nr. of sites: ",Nsites)
    println(f)
    close(f)
end

function write_output(path,g,observable)
    text=' '
    for b=1:length(observable)
        text*=string(observable[b],' ')
    end
    f=open(path,"a")
    println(f,round(g,digits=4)," ",text)
    close(f)
end

function write_text(path,g,txt)
    f=open(path,"a")
    println(f,round(g,digits=4)," ",txt)
    close(f)
end

###
###save data from Observers, code by Wladislaw Krinitsin
# util function to save data in obs
function savedata(name::String, obs; kwargs...)
  name == "" && return
  h5open(name*".h5", "w") do file    
    for pair in kwargs
      file["variables/"*String(pair[1])] = pair[2]
    end

    # iterate through the fields of obs and append the data to the dataframe
    for n in names(obs)
        # println(file)
        if obs[1,n] isa Array # if vector of vectors, reduce to matrix
            tmp = mapreduce(permutedims, vcat, obs[!,n])
            # println(tmp)
            write_dataset(file, "data/"*n, tmp)
        else
            write_dataset(file, "data/"*n, obs[!,n])
        end
    end
  end
end
  
# util function to save data in obs
function saveparams(name::String, params)
  name == "" && return
  h5open(name*".h5", "cw") do file    
    create_group(file, "params")
    for (name,val) in params
      file["params"][name] = val
    end 
  end
end# util function to save data in obs


function MPS_to_ITensorNetwork(mps)
    return ITensorNetwork([mps[v] for v in eachindex(mps)])
end
function TTN_to_MPS(tn)
    return MPS(collect(ITensorNetworks.vertex_data(tn)))
    # TTN([tmp[v] for v in eachindex(tmp)])
end
#end