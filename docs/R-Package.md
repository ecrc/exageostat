**ExaGeoStatR**

exageostatR is an R package for ExaGeoStat framework: a parallel high performance unified framework for geostatistics on manycore systems

**Getting Started**

**install dependencies**

1. Portable Hardware Locality (hwloc).
2. NLopt.
3. GNU Scientific Library (GSL).
4. The StarPU runtime system.
5. The CHAMELEON dense linear algebra software.
6. Easy installation of the above packages is available by using install_cpu.sh

git clone exageostatR repo

git clone https://github.com/ecrc/exageostatR.git

**Update submodules**

git submodule init

git submodule update

**Build R package**

mv ./r-wrappers ./src

mv ./examples ./src

mv ./src/Makefile-shlib ./src/Makefile

R CMD check exageostatR

R CMD BUILD exageostatR

**Use ExaGeoStatR**

install.packages(repos=NULL, "exageostatR_0.1.0.tar.gz")

library("exageostatR")

Possibilities of ExaGeoStat

**Operations:**

1. Generate synthetic spatial datasets (i.e., locations & environmental measurements).
2. Maximum likelihood estimation using dense covariance matrices. 
3. Task-based programmin model and dynamic runtime system using CHAMELEON & StarPU

**Tutorial**

A more detailed description can be accessible here

For example, to search for data scientist jobs in London:

\code{.R}
library("exageostat")
library("RhpcBLASctl")
#Inputs
theta1 = 1       # initial variance
theta2 = 0.1     # initial smoothness
theta3 = 0.5     # initial range
computation = 0  # exact computation
dmetric = 0      # 'ed'  Euclidean distance
n=1600           # n*n locations grid 
gpus=0           # number of underlying GPUs
ts=320           # tile_size:  change it could improve the performance. No fixed value can be given
p_grid=1         # more than 1 in the case of distributed systems 
q_grid=1         # more than 1 in the case of distributed systems ( usually equals to p_grid)
clb = vector(mode="numeric",length = 3)    #optimization lower bounds
cub = vector(mode="numeric",length = 3)    #optimization upper bounds
theta_out = vector(mode="numeric",length = 3)    # parameter vector output
clb=as.numeric(c("0.01","0.01","0.01"))
globalveclen =  3*n
cub=as.numeric(c("5","5","5"))
vecs_out = vector(mode="numeric",length = globalveclen)     #Z measurements of  locations
vecs_out[1:globalveclen] = -1.99
theta_out[1:3]= -1.99
#Generate Z observation vector
vecs_out = rexageostat_gen_zR(n, ncores, gpus, ts, p_grid, q_grid, theta1, theta2, theta3, computation, dmetric, globalveclen)
#Estimate MLE parameters
theta_out = rexageostat_likelihoodR(n, ncores, gpus, ts, p_grid, q_grid,  vecs_out[1:n],  vecs_out[n+1:(2*n)],  vecs_out[(2*n+1):(3*n)], clb, cub, computation, dmetric)
  \endcode

