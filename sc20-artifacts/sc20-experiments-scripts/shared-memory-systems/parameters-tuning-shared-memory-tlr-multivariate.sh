#!/bin/bash

#dense

ncores=55 #intel Skylake example

for n in     28900  40000  55225 59536  
do
                        for k in 1 2  3 
                        do
                                for dts in 320 460 560 760 960
                                do
                                        ../build/examples/synthetic_dmle_test --test          --N=$n          --dts=$dts          --ncores=$ncores          --computation=exact          --kernel=1:1:0.03:0.5:1:0.9          --ikernel=1:1:0.03:0.5:1:0.9          --olb=0.01:0.01:0.01:0.01:0.01:0.01          --oub=2:2:2:2:2:1          --zvecs=1       --kernel_fun=bivariate_matern_parsimonious     
                                done
                        done
done



#LR5
for n in     28900  40000  55225 59536
do
	for rank in  275 300 325 350 
	do
		for lts in 960 1000 1200 1400 1600 
		do
			for k in 1 23
			do
        			for dts in 560   
        			do
					../build/examples/synthetic_dmle_test --test          --N=$n          --dts=$dts          --ncores=$ncores          --computation=lr_approx          --kernel=1:1:0.03:0.5:1:0.9          --ikernel=1:1:0.03:0.5:1:0.9          --olb=0.01:0.01:0.01:0.01:0.01:0.01          --oub=2:2:2:2:2:1          --zvecs=1       --kernel_fun=bivariate_matern_parsimonious    --acc=5  --maxrank=$rank --lts=$lts
        			done
			done
		done
	done	
done


#LR7
for n in     28900  40000  55225 59536
do
       for rank in  325 350  375 400 425 450 475 500
       do
               for lts in 960 1000 1200 1400 1600 1800 2000 2400
               do
                       for k in 1 2 3
                       do
                               for dts in 560
                               do
                                       ../build/examples/synthetic_dmle_test --test          --N=$n          --dts=$dts          --ncores=$ncores          --computation=lr_approx          --kernel=1:1:0.03:0.5:1:0.9          --ikernel=1:1:0.03:0.5:1:0.9          --olb=0.01:0.01:0.01:0.01:0.01:0.01          --oub=2:2:2:2:2:1          --zvecs=1       --kernel_fun=bivariate_matern_parsimonious    --acc=7  --maxrank=$rank --lts=$lts
                               done
                       done
               done
       done
done



#lR9
for n in     28900  40000  55225 59536
do
       for rank in  600 625 650 675 700 
       do
              for lts in 960 1000 1200 1400 1600 1800 2000 2400
               do
                       for k in 1 2 3
                       do
                               for dts in 560
                               do
                                       ../build/examples/synthetic_dmle_test --test          --N=$n          --dts=$dts         --ncores=$ncores          --computation=lr_approx          --kernel=1:1:0.03:0.5:1:0.9          --ikernel=1:1:0.03:0.5:1:0.9          --olb=0.01:0.01:0.01:0.01:0.01:0.01          --oub=2:2:2:2:2:1          --zvecs=1       --kernel_fun=bivariate_matern_parsimonious    --acc=9  --maxrank=$rank --lts=$lts
                               done
                       done
               done
       done
done
