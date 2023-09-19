mutable struct DemoObserver <: AbstractObserver
        energy_tol::Float64
        last_energy::Float64

        DemoObserver(energy_tol=0.0) = new(energy_tol,1000.0)
end

function ITensors.checkdone!(o::DemoObserver;kwargs...)
        sw = kwargs[:sweep]
        energy = kwargs[:energy]
        wf=kwargs[:psi]
        maxbond=maxlinkdim(wf)
        en_error=abs(energy-o.last_energy)/abs(energy)

        fle=open(joinpath(output_path, "log"),"a")
        @printf(fle,"Sweep %i E= %.10e error= %.10e MaxBondDim= %i \n",sw,energy,en_error,maxbond)
        if en_error < o.energy_tol
                println(fle,"Convergence reached after sweep $sw")
                #println(fle,"Max. bond dimension: ",maxbond)
                #println(fle,"E(final)= ",energy)
                close(fle)
                return true
        end
        close(fle)

        # Otherwise, update last_energy and keep going
        o.last_energy = energy
        return false
end
