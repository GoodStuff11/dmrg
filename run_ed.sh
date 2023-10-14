#!/bin/bash
#SBATCH --account=def-pnroy 
#SBATCH --mem=80G
#SBATCH --time=2-00:00

module load julia/1.9.1

julia --project=. ed/ed.jl
