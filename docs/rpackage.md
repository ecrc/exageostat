**ExaGeoStatR 1.0.0**

exageostatR is an R package for ExaGeoStat framework: a parallel high performance unified framework for geostatistics on manycore systems

**Getting Started**

**install dependencies**

1. Portable Hardware Locality (hwloc).
2. NLopt.
3. GNU Scientific Library (GSL).
4. The StarPU runtime system.
5. The CHAMELEON dense linear algebra software.
6. STARS-H, a high performance parallel open-source package of Software for 
   Testing Accuracy, Reliability and Scalability of Hierarchical computations
7. The Hierarchical Computations on Manycore Architectures (HiCMA) library.


git clone exageostatR repo

git clone https://github.com/ecrc/exageostatR.git

**Update submodules**

git submodule init

git submodule update

**Build R package**


R CMD check exageostatR

R CMD BUILD exageostatR

**Use ExaGeoStatR**

install.packages(repos=NULL, "exageostatR_1.0.0.tar.gz")

library("exageostat")

Possibilities of ExaGeoStat

**Operations:**

1. Generate synthetic spatial datasets (i.e., locations & environmental measurements).
2. Maximum likelihood estimation using dense, Diagonal Super Tile (DST), and Tile Low-Rank (TLR)  covariance matrices. 
3. Task-based programmin model and dynamic runtime system using HiCMA, CHAMELEON, and StarPU.

**Installation User Guide and Examples**
A complete Installation User Guide can be accessiable [here](https://github.com/ecrc/exageostatr/exageostatr.pdf).
A more detailed description and examples be accessible [here](https://github.com/ecrc/exageostatr).


