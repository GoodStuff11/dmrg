#!/bin/bash
#SBATCH --account=rrg-pnroy 
#SBATCH --mem=10G
#SBATCH --time=2-00:00

module load julia/1.9.3
cd dmrg

# if [ $# -eq 0 ] ; then
#     julia --project=.. ./dmrg_run.jl  
# else
#     julia --project=.. ./dmrg_run.jl gstart $1 ParitySymmetry $2
# fi

julia --project=../ ./dmrg_analysis.jl  