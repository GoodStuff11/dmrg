for N in 5 10 15 20 50 80 100; do
    for cutoff in 1e-6 1e-7 1e-8; do 
        sbatch run_tdvp.sh $N $cutoff
    done
done
