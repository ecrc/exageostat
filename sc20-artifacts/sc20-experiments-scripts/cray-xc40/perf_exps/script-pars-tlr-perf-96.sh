#!/bin/bash
#SBATCH --job-name=exact_%j
#SBATCH --output=./output_scala/exact_96_%A.out
#SBATCH -A k1339
#SBATCH --nodes=96
#SBATCH --ntasks=96
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=08:30:00

export STARPU_SCHED=eager

export STARPU_WATCHDOG_TIMEOUT=1000000000
export STARPU_WATCHDOG_CRASH=1


#x=10000
#y=9000
#export STARPU_LIMIT_MAX_SUBMITTED_TASKS=$x
#export STARPU_LIMIT_MIN_SUBMITTED_TASKS=$y

n=$1
dts=$2



echo "srun ./build/examples/synthetic_dmle_test --test          --N=$n          --dts=$dts          --ncores=30          --computation=exact          --kernel=1:1:0.03:0.5:1:0.9          --ikernel=1:1:0.03:0.5:1:0.9          --olb=0.01:0.01:0.01:0.01:0.01:0.01          --oub=2:2:2:2:2:1          --zvecs=1       --kernel_fun=bivariate_matern_parsimonious      --p=12 --q=8 --opt_iter=1"


 srun ./build/examples/synthetic_dmle_test --test          --N=$n          --dts=$dts          --ncores=30          --computation=exact          --kernel=1:1:0.03:0.5:1:0.9          --ikernel=1:1:0.03:0.5:1:0.9          --olb=0.01:0.01:0.01:0.01:0.01:0.01          --oub=2:2:2:2:2:1          --zvecs=1       --kernel_fun=bivariate_matern_parsimonious      --p=12  --q=8 --opt_iter=1

