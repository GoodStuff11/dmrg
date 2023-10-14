#!/bin/bash
#SBATCH --account=def-pnroy 
#SBATCH --mem=80000M
#SBATCH --time=2-00:00

module load julia/1.9.1
cd dmrg_tdvp

if [ $# -eq 0 ] ; then
    julia --project=../ ./dmrg.jl 
else
    julia --project=../ ./dmrg.jl $1 $2 $3 $4 $5 $6
fi