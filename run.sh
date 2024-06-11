#!/bin/bash
#SBATCH --account=rrg-pnroy 
#SBATCH --mem=32G
#SBATCH --time=2-00:00

module load julia/1.9.1
cd dmrg_tdvp

if [ $# -eq 0 ] ; then
    julia --project=../ ./tdvp.jl 
else
    julia --project=../ ./tdvp.jl $1
fi