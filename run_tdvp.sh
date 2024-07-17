#!/bin/bash
#SBATCH --account=rrg-pnroy 
#SBATCH --mem=60G
#SBATCH --time=2-00:00

module load julia/1.9.3
cd dmrg

if [ $# -eq 0 ] ; then
    julia --project=.. ./tdvp.jl  
else
    julia --project=.. ./tdvp.jl Nsites $1 ecutoff $2
fi
