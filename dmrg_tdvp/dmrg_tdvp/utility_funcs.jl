module utility_funcs

using Printf

export create_file, write_output, write_text

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

end