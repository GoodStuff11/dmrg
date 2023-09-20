module utility_funcs

using Printf
export create_file, write_output, write_text, Logger

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

mutable struct Logger
    path::String
    filename::String

    function Logger(path::String)
        new(path, "log")
    end
    function Logger(path::String, filename::String)
        new(path, filename)
    end
end

function Base.:<(log::Logger, strings::Vector{Any})
    file = open(joinpath(log.path, log.filename),"a")
    for s in strings
        if s isa Tuple
            println(file, s...)
        else
            println(file,s)
        end
    end
    close(file)
end

function Base.:<(log::Logger, string::Union{Any, Tuple{Any}})
    file = open(joinpath(log.path, log.filename),"a")
    if string isa Tuple
        println(file,string...)
    else
        println(file,string)
    end
    close(file)
end


end

