###################################################################################### Real Dataset
#ncores=4;

#for ((ts=300;ts<=900;ts+=100))
#do
#     for ((b=ncores; b<=32; b*=2))
#     do

#         for((c=1; c<=1; c+=1))
#         do
#		./MLE --ncores=$b  --ts=$ts  --locs_file=./NCAR_pinfill_temperatures/METAinfo --obs_dir=./NCAR_pinfill_temperatures --timestamp=$c --async > ./results/$a-$ts-$b-atp-matern-NCAR_pinfill_temperatures.txt
 #               ./MLE --ncores=$b  --ts=$ts  --locs_file=./NCAR_pinfill_temperatures/METAinfo --obs_dir=./NCAR_pinfill_temperatures --timestamp=$c --chameleon > ./results/$a-$ts-$b-stc-matern-NCAR_pinfill_temperatures.txt
  #              ./MLE --ncores=$b  --ts=$ts  --locs_file=./NCAR_pinfill_precipitation/METAinfo --obs_dir=./NCAR_pinfill_precipitation --timestamp=$c --async > ./results/$a-$ts-$b-atp-matern-NCAR_pinfill_precipitation.txt
   #             ./MLE --ncores=$b  --ts=$ts  --locs_file=./NCAR_pinfill_precipitation/METAinfo --obs_dir=./NCAR_pinfill_precipitation --timestamp=$c --chameleon > ./results/$a-$ts-$b-stc-matern-NCAR_pinfill_precipitation.txt


#	done
 #    done
#done     


###############################################################################################################################################################################################################
#N=14000;

#for ((a=N; a <= 46000 ; a+=8000))  # Double parentheses, and "LIMIT" with no "$".
#do
#		last_ts=$a/28;
#		for ((ts=200;ts<=$last_ts;ts+=100))
#		do
#			if [ $a -le 30000 ]
  #                      then	
 #                          ts=300;
#			   last_ts=200;		
 #                       fi
			#./MLE --test --N=$a --ts=$ts --ncores=$b  --async  --kernel=1:?:0.5 --ikernel=1:0.1:0.5 > ./results/$a-$ts-$b-atp-exp.txt
        #	       numactl --interleave=all ./MLE --test --N=$a --ts=$ts --ncores=28  --async --kernel=1:?:0.5 --ikernel=1:0.1:0.5 > ./results/$a-$ts-$b-stc-exp.txt
	               # ./MLE --test --N=$a --ts=$ts --ncores=$b  --async  --kernel=?:?:? --ikernel=1:0.05:1 > ./results/$a-$ts-$b-atp-matern.txt
        	       # ./MLE --test --N=$a --ts=$ts --ncores=$b  --chameleon  --kernel=?:?:? --ikernel=1:0.05:1 > ./results/$a-$ts-$b-stc-matern.txt

#		done
   
#done



###############################################################################################################################################################################################################
N=100000;
ncores1=26;

for ((a=N; a <= 150000 ; a+=25000))  # Double parentheses, and "LIMIT" with no "$".
do
echo N=$a
echo ncores=$ncores1
                for ((ts=400; ts<=(($a/$ncores1)); ts+=50))
                do
echo ts=$ts

              STARPU_SILENT=1  numactl --interleave=all ./MLE --test --N=$a --ts=$ts --ncores=$ncores1  --chameleon  --kernel=1:?:0.5 --ikernel=1:0.1:0.5 > ./results/$a-$ts-$b-atp-exp.txt
              #  ./MLE --test --N=$a --ts=$ts --ncores=$b  --chameleon  --kernel=1:?:0.5 --ikernel=1:0.1:0.5 > ./results/$a-$ts-$b-stc-exp.txt
               # ./MLE --test --N=$a --ts=$ts --ncores=$b  --async  --kernel=?:?:? --ikernel=1:0.05:1 > ./results/$a-$ts-$b-atp-matern.txt
               # ./MLE --test --N=$a --ts=$ts --ncores=$b  --chameleon  --kernel=?:?:? --ikernel=1:0.05:1 > ./results/$a-$ts-$b-stc-matern.txt
                done
   
done



###############################################################################################################################################################################################################
#N=30000;

#for ((a=N; a <= 120000 ; a+=30000))  # Double parentheses, and "LIMIT" with no "$".
#do

#last_ts=($a/64)
 #               for ((ts=300;ts<=last_ts;ts+=100))
  #              do
#	if [ $a -le 14000 ]
 #                       then	
 #                          ts=300;
#			   last_ts=200;		
 #                       fi
#	if [ $a -le 22000 ]
 #                       then	
  #                         ts=500;
#			   last_ts=200;		
 #                       fi
#	if [ $a -le 30000 ]
 #                       then	
  #                         ts=700;
#			   last_ts=200;		
 #                       fi
#	if [ $a -le 38000 ]
 #                       then	
  #                         ts=700;
#			   last_ts=200;		
 #                       fi
#	if [ $a -le 46000 ]
 #                       then	
  #                         ts=700;
#			   last_ts=200;		
 #                       fi

                #./MLE --test --N=$a --ts=$ts --ncores=$b  --async  --kernel=1:?:0.5 --ikernel=1:0.1:0.5 > ./results/$a-$ts-$b-atp-exp.txt
#                numactl --interleave=all    ./MLE --test --N=$a --ts=$ts --ncores=64  --chameleon  --kernel=?:?:? --ikernel=1:0.05:1 > ./results/$a-$ts-$b-stc-exp.txt
               # ./MLE --test --N=$a --ts=$ts --ncores=$b  --async  --kernel=?:?:? --ikernel=1:0.05:1 > ./results/$a-$ts-$b-atp-matern.txt
               # ./MLE --test --N=$a --ts=$ts --ncores=$b  --chameleon  --kernel=?:?:? --ikernel=1:0.05:1 > ./results/$a-$ts-$b-stc-matern.txt

                done
         
done



###############################################################################################################################################################################################################
#N=14000;

##for ((a=N; a <= 46000 ; a+=8000))  # Double parentheses, and "LIMIT" with no "$".
#do

 #             last_ts=($a/28)
  #              for ((ts=200;ts<=last_ts;ts+=100))
   #             do
#	if [ $a -le 14000 ]
 #                       then	
  #                         ts=300;
#			   last_ts=200;		
 #                       fi
#	if [ $a -le 22000 ]
 #                       then	
  #                         ts=500;
#			   last_ts=200;		
 #                       fi
#	if [ $a -le 30000 ]
 #                       then	
  #                         ts=700;
#			   last_ts=200;		
 #                       fi
#	if [ $a -le 38000 ]
 #                       then	
  #                         ts=200;
#			   last_ts=100;		
 #                       fi
#	if [ $a -le 46000 ]
 #                       then	
  #                         ts=200;
#			   last_ts=100;		
 #                       fi
#


 #             numactl --interleave=all    ./MLE --test --N=$a --ts=$ts --ncores=36  --async  --kernel=?:?:? --ikernel=1:0.05:1 > ./results/$a-$ts-$b-atp-exp.txt
              #  ./MLE --test --N=$a --ts=$ts --ncores=$b  --chameleon  --kernel=1:?:0.5 --ikernel=1:0.1:0.5 > ./results/$a-$ts-$b-stc-exp.txt
               # ./MLE --test --N=$a --ts=$ts --ncores=$b  --async  --kernel=?:?:? --ikernel=1:0.05:1 > ./results/$a-$ts-$b-atp-matern.txt
               # ./MLE --test --N=$a --ts=$ts --ncores=$b  --chameleon  --kernel=?:?:? --ikernel=1:0.05:1 > ./results/$a-$ts-$b-stc-matern.txt

  #              done

#done
