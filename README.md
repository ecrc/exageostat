What is ExaGeoStat?
================

The **Exascale GeoStatistics** project (ExaGeoStat) is a parallel high performance unified framework for computational geostatistics on many-core systems. The project aims at optimizing the likelihood function for a given spatial data to provide an efficient way to predict missing observations in the context of climate/weather forecasting applications. This machine learning framework proposes a unified simulation code structure to target various hardware architectures, from commodity x86 to GPU accelerator-based shared and distributed-memory systems. ExaGeoStat enables statisticians to tackle computationally challenging scientific problems at large-scale, while abstracting the hardware complexity through state-of-the-art high performance linear algebra software libraries.



Vision of ExaGeoStat
=================

ExaGeoStat is a collaboration between the KAUST Spatial Statistics group and the Extreme Computing Research
Center (ECRC). Its contribution lies not in a new algorithm nor in a new data set,
but in demonstrating the routine use of the larger data sets becoming available to geospatial
statisticians, thanks to the implementation of state-of-the-art statistical algorithms on
high-performance computing (HPC) hardware. 

We have built a standalone software framework (ExaGeoStat) that is able to run on a variety
of hardware resources, including GPUs and massively distributed systems such as Shaheen-II,
KAUST's Cray XC40 supercomputer, ORNL Summit (OLCF-4) supercomputer, and Riken Fugaku supercomputer, 
to create a statistical model to predict environmental data (i.e., temperature, flow rates, soil moisture,
 wind speed, air pollution, etc.) at spatial locations on which data
is missing, and to exploit large amounts of data to reduce the effect of individual measurement
errors.  The best-known methods for such statistical processing have a cost that grows rapidly
in the size of the data set, namely, in proportion to its cube, or third power.  Thus, increasing
the size of the data set by a factor ten drives up the cost of the computation by a factor of
a thousand, while simultaneously driving up the memory requirements by a factor a hundred.  

For instance, according to this cubic growth in complexity, a computation that requires one
minute would require nearly 17 hours on a data set just ten times larger. This creates a
computational strain on standard statistics software, for which contemporary data sizes
were not anticipated; and even if possible, it puts the computation beyond the interactive
attention span of the analyst. Parallelism (assigning thousands of processors to the
single task) and Moore's Law allow leading-edge computers to handle such "big data" 
with ease, but the software bridge must be built.  Furthermore, the software interface
must resemble the interactive one with which working statisticians are familiar.

To summarize, the combination of emerging computing capabilities and emerging datasets
promises significant advances in statistical analyses of environmental and many other
phenomena.  Such cross-disciplinary advances are natural at KAUST, which is why this
relatively low-hanging fruit was ours to harvest earliest. Our roadmap takes now ExaGeoStat 
a step further on the algorithmic side by integrating tile low-rank matrix approximation.
This low-rank matrix approximation permits to exploit the data sparisty of the operator
with a user-controlled numerical accuracy. This further expands practical problem sizes for
statisticians with modest computational resources.


Current Version: 1.2.1
======================
Current Features of ExaGeoStat
======================
Operations:
1. Large-scale synthetic matrix generation in dense.
2. Univariate/bivariate/trivariate Gaussian space modeling and prediction using dense,
 Tile Low-Rank (TLR), Diagonal Super-Tile (DST) approximations.
3. Univariate/bivariate Gaussian space/time modeling and prediction using dense,
 Tile Low-Rank (TLR), Diagonal Super-Tile (DST) approximations.
4. Univariate Gaussian and Tukey g-and-h non-Gaussian space modeling using dense,
 Tile Low-Rank (TLR), Diagonal Super-Tile (DST) approximations.
5. Univariate and bivariate parameter estimation assessment using MLOE/MMOM criteria.

Supported Covariance Functions:
======================
1. Univariate Matérn (Gaussian/Stationary)
2. Univariate Matérn with Nugget (Gaussian/Stationary)
3. Flexible Bivariate Matérn (Gaussian/Stationary)
4. Parsimonious Bivariate Matérn (Gaussian/Stationary)
5. Parsimonious trivariate Matérn (Gaussian/Stationary)
6. Univariate Space/Time Matérn (Gaussian/Stationary)
7. Bivariate Space/Time Matérn (Gaussian/Stationary)
8. Tukey g-and-h Univariate Matérn (non-Gaussian/Stationary)
9. Tukey g-and-h Univariate Power Exponential (non-Gaussian/Stationary)

Programming models:
1.  MPI
2.  Task-based programming models

External libraries:
1.  NLOPT
2.  GSL
3.  HWLOC
4.  StarPU dynamic runtime system 
5.  HiCMA
6.  Stars-H
7.  Chameleon

Installation
============

Installation requires at least **CMake** of version 2.8.12. To build ExaGeoStat,
please follow these instructions:

1.  Get  from git repository

        git clone git@github.com:ecrc/exageostat

    or

        git clone https://github.com/ecrc/exageostat

2.  Go into ExaGeoStat folder

        cd exageostat

3.  Get submodules

        git submodule update --init --recursive

4.  Create build directory and go there

        mkdir -p build && cd build

5.  Use CMake to get all the dependencies

        cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/  -DExaGeoStat_SCHED_STARPU=ON   -DEXAGEOSTAT_USE_NETCDF=ON -DEXAGEOSTAT_USE_HICMA=ON -DEXAGEOSTAT_USE_CHAMELEON=ON -DEXAGEOSTAT_INSTALL_DEPS=ON -DBUILD_SHARED_LIBS=ON

6.  Build EXAGEOSTAT

        make -j

7.  Build local documentation (optional)

        make docs

8.  Install EXAGEOSTAT

        make install

9. Add line

        export PKG_CONFIG_PATH=/path/to/install/lib/pkgconfig:$PKG_CONFIG_PATH

    to your .bashrc file.

Now you can use pkg-config executable to collect compiler and linker flags for
EXAGEOSTAT.

Check docs/install.md for more installations details.

References
==========
1. Abdulah, Sameh, Hatem Ltaief, Ying Sun, Marc G. Genton, and David E. Keyes. "ExaGeoStat: A high performance unified software for geostatistics on manycore systems." IEEE Transactions on Parallel and Distributed Systems 29, no. 12 (2018): 2771-2784.

2. Abdulah, Sameh, Hatem Ltaief, Ying Sun, Marc G. Genton, and David E. Keyes. "Parallel approximation of the maximum likelihood estimation for the prediction of large-scale geostatistics simulations." In 2018 IEEE international conference on cluster computing (CLUSTER), pp. 98-108. IEEE, 2018.

3. Abdulah, Sameh, Hatem Ltaief, Ying Sun, Marc G. Genton, and David E. Keyes. "Geostatistical modeling and prediction using mixed precision tile Cholesky factorization." In 2019 IEEE 26th international conference on high performance computing, data, and analytics (HiPC), pp. 152-162. IEEE, 2019.

4. Salvana, Mary Lai O., Sameh Abdulah, Huang Huang, Hatem Ltaief, Ying Sun, Marc G. Genton, and David E. Keyes. "High performance multivariate geospatial statistics on manycore systems." IEEE Transactions on Parallel and Distributed Systems 32, no. 11 (2021): 2719-2733.

5. Salvaña, Mary Lai O., Sameh Abdulah, Hatem Ltaief, Ying Sun, Marc G. Genton, and David E. Keyes. "Parallel Space-Time Likelihood Optimization for Air Pollution Prediction on Large-Scale Systems." In Proceedings of the Platform for Advanced Scientific Computing Conference (PASC '22). Association for Computing Machinery, New York, NY, USA, Article 17, 1–11. ACM, 2022.

6. Abdulah, Sameh, Qinglei Cao, Yu Pei, George Bosilca, Jack Dongarra, Marc G. Genton, David E. Keyes, Hatem Ltaief, and Ying Sun. "Accelerating geostatistical modeling and prediction with mixed-precision computations: A high-productivity approach with PaRSEC." IEEE Transactions on Parallel and Distributed Systems 33, no. 4 (2021): 964-976.

7. Mondal, Sagnik, Sameh Abdulah, Hatem Ltaief, Ying Sun, Marc G. Genton, and David E. Keyes. "Parallel Approximations of the Tukey g-and-h Likelihoods and Predictions for Non-Gaussian Geostatistics." 2022 IEEE International Parallel and Distributed Processing Symposium (IPDPS), Lyon, France, 2022, pp. 379-389. IEEE, 2022.

8. Cao, Qinglei, Sameh Abdulah, Rabab Alomairy, Yu Pei, Pratik Nag, George Bosilca, Jack Dongarra et al. "Reshaping geostatistical modeling and prediction for extreme-scale environmental applications." In 2022 SC22: International Conference for High Performance Computing, Networking, Storage and Analysis (SC), pp. 13-24. IEEE Computer Society, 2022. (ACM GORDON BELL PRIZE Finalist).

Handout
=======

More information can be found in this handout:
![Handout](docs/handout.jpg)
