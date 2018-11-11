What is ExaGeoStat?
================

The **Exascale GeoStatistics** project (ExaGeoStat) is a parallel high performance unified framework for computational geostatistics on many-core systems. The project aims at optimizing the likelihood function for a given spatial data to provide an efficient way to predict missing observations in the context of climate/weather forecasting applications. This machine learning framework proposes a unified simulation code structure to target various hardware architectures, from commodity x86 to GPU accelerator-based shared and distributed-memory systems. ExaGeoStat enables statisticians to tackle computationally challenging scientific problems at large-scale, while abstracting the hardware complexity, through state-of-the-art high performance linear algebra software libraries.



Vision of ExaGeoStat
=================

ExaGeoStat is a collaboration between the KAUST Statistics group and the Extreme Computing Research
Center (ECRC). Its contribution lies not in a new algorithm nor in a new data set,
but in demonstrating the routine use of the larger data sets becoming available to geospatial
statisticians, thanks to the implementation of state-of-the-art statistical algorithms on
high-performance computing (HPC) hardware. 

We have built a standalone software framework (ExaGeoStat) that is able to run on a variety
of hardware resources, including GPUs and massive distributed systems such as Shaheen,
KAUST's Cray XC40 supercomputer, and to create a statistical model to predict environmental data
(i.e., temperature, flow rates, soil moisture, wind speed, etc.) at spatial locations on which data
is missing, and to exploit large amounts of data to reduce the effect of individual measurement
errors.  The best-known methods for such statistical processing have a cost that grows rapidly
in the size of the data set, namely, in proportion to its cube, or third power.  Thus, increasing
the size of data set by a factor ten drives up the cost of the computation by a factor of
a thousand, while simultaneously driving up the memory requirements by a factor o hundred.  

For instance, according to this cubic growth in complexity, a computation that requires one
minute would require nearly 17 hours on a data set just ten times larger. This creates a
computational strain on standard statistics software, for which contemporary data sizes
were not anticipated; and even if possible, it puts the computation beyond the interactive
attention span of the analyst. Parallelism (assigning thousands of processors to the
single task) and Moore's Law allow leading edge computers to handle such "big data" 
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


Current Version: 1.0.0
======================

Current Features of ExaGeoStat
======================
Operations:
1.  Large-scale synthetic matrix generation.
2.  Maximum likelihood estimation from dense and Tile Low-Rank (TLR) covariance matrices.
3.  Predicting large-scale unknown measures in predefined geo-spatial locations.

Programming models:
1.  MPI
2.  Task-based programming models

External libraries:
1.  StarPU dynamic runtime system 
2.  HiCMA
3.  Chameleon

Installation
============

Installation requires at least **CMake** of version 3.2.3. To build ExaGeoStat,
please follow these instructions:

1.  Get  from git repository

        git clone git@github.com:ecrc/exageostat

    or

        git clone https://github.com/ecrc/exageostat

2.  Go into exageostat folder

        cd exageostat

3.  Get submodules

        git submodule update --init --recursive

4.  Create build directory and go there

        mkdir build && cd build

5.  Use CMake to get all the dependencies

        cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/  -DEXAGEOSTAT_SCHED_STARPU=ON   -DEXAGEOSTAT_USE_NETCDF=ON -DEXAGEOSTAT_USE_HICMA=ON


6.  Build EXAGEOSTAT

        make -j

8.  Build local documentation (optional)

        make docs

9.  Install EXAGEOSTAT

        make install

10. Add line

        export PKG_CONFIG_PATH=/path/to/install/lib/pkgconfig:$PKG_CONFIG_PATH

    to your .bashrc file.

Now you can use pkg-config executable to collect compiler and linker flags for
EXAGEOSTAT.


References
==========
1. Sameh Abdulah, Hatem Ltaief, Ying Sun, Marc G. Genton, and David Keyes. "ExaGeoStat: A High Performance Unified Software for Geostatistics on Manycore Systems," IEEE Transactions on Parallel and Distributed Systems (2018).

2. Sameh Abdulah, Hatem Ltaief, Ying Sun, Marc G. Genton, and David Keyes. "Parallel Approximation of the Maximum Likelihood Estimation for the Prediction of Large-Scale Geostatistics Simulations," IEEE Cluster Conference, Belfast, UK, Septemeber, 2018.


Handout
=======

More information can be found in this handout:
![Handout](docs/handout.jpg)
