#!/bin/bash
#SBATCH --account=rrg-pnroy 
#SBATCH --mem=32G
#SBATCH --time=2-00:00

module load julia/1.9.3

julia --project=. ed/ed.jl
