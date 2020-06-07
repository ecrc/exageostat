#!/bin/bash
##dense

ncores=55 # Intel Skylake example

for n in 22500 28900 34225 40000 46225 51076   55225 63001 
do
	for k in 1 2 3
	do
		predict_values=100
		echo $predict_values
        	for dts in 560   
        	do
		../build/examples/synthetic_dmle_test --test          --N=$n          --dts=$dts          --ncores=$ncores          --computation=exact          --kernel=?:?:?:?:?:? --opt_iter=3     --ikernel=1:1:0.1:0.08:0.5:0.8         --olb=0.01:0.01:0.01:0.01:0.01          --oub=5:5:5:5:5:5          --zvecs=1        --kernel_fun=bivariate_matern_parsimonious   --predict=$predict_values  --mloe_mmom
        	done
	done	
done

