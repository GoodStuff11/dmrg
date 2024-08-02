# k=0
# for i in $(seq 0.1 0.1 3); do
#     ((k++))
#     sbatch run.sh "data_$k" $i 0.1 1 5 6
# done


for Nspec in 9 11 13; do
    for g in 0.25 0.5; do
        # echo $Nsites $mmax
        sbatch run_dmrg.sh $g even even $Nspec 100
        sbatch run_dmrg.sh $g odd odd $Nspec 100
        sbatch run_dmrg.sh $g even odd $Nspec 100
        sbatch run_dmrg.sh $g odd even $Nspec 100
    done
done