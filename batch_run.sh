k=0
for i in $(seq 0.1 0.1 3); do
    ((k++))
    sbatch run.sh "data_$k" $i 0 1 5
done