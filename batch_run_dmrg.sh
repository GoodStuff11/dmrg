# k=0
# for i in $(seq 0.1 0.1 3); do
#     ((k++))
#     sbatch run.sh "data_$k" $i 0.1 1 5 6
# done

k=0

for g in $(seq 0.1 0.4 1.9); do
    ((k++))
    # echo $Nsites $mmax
    sbatch run_dmrg.sh $g even even
    sbatch run_dmrg.sh $g odd odd
    sbatch run_dmrg.sh $g even odd
    sbatch run_dmrg.sh $g odd even
done
