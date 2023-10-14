# k=0
# for i in $(seq 0.1 0.1 3); do
#     ((k++))
#     sbatch run.sh "data_$k" $i 0.1 1 5 6
# done

k=0
for Nsites in $(seq 4 1 12); do
    for mmax in $(seq 1 1 5); do
        ((k++))
        # echo $Nsites $mmax
        sbatch run.sh "data_$k" 0.1 0.1 30 $mmax $Nsites
    done
done