#!/bin/bash
export OMP_NUM_THREADS=36
#for ((W=26624; W <= 210000 ; W+=8192))  # Double parentheses, and "LIMIT" with no "$".
#	do	
#	numactl --interleave=all	./LAPACKE_MLE 28 $W 0 1; 
#done




numactl --interleave=all        ./LAPACKE_MLE 28 1600 0 1 1 0.1 0.5;
numactl --interleave=all        ./LAPACKE_MLE 28 51076 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 55225 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 59049 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 63001 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 67081 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 71289 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 75625 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 79524 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 83521 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 87616 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 91809 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 96100 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 99856 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 104329 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 108241 0 1 1 0.1 0.5
numactl --interleave=all        ./LAPACKE_MLE 28 112225 0 1 1 0.1 0.5


