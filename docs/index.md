What is ExaGeoStat? {#mainpage}
===================

The **Exascale Geostatistics** project (ExaGeoStat) is a parallel high performance 
unified framework for geostatistics on manycore systems. The framework aims at optimizing the
likelihood function for a given spatial data in order to provide an efficient methodology 
for predicting missing observations, e.g., in the context of climate and weather forecasting. 
The framework targets many-core shared and distributed-memory systems, composed of x86 
and GPUs.

### Vision of ExaGeoStat ###

ExaGeoStat is a collaboration between the KAUST Statistics group and the Extreme Computing Research
Center (ECRC). Its contribution lies not in a new algorithm nor in a new data set,
but in demonstrating routine use of the larger data sets becoming available to geospatial
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
relatively low-hanging fruit was ours to harvest earliest. Our roadmap is take ExaGeoStat 
a step further on the algorithmic side by introducing a new type of approximation
(hierarchical tile low-rank approximation) that is systematic and whose errors can
be understood and controlled, in order to further expand practical problem sizes for
even modest computational resources, since most working statisticians lack access to
systems like Shaheen.

### Current Features of ExaGeoStat ###

Operations:
1.  Generating synthetic datasets that mimic climate/weather model (i.e., generate longitudes, latitudes and observations),
2.  Evaluating Maximum Likelihood function for both real and synthetic datasets using exact, Diagonal Super Tile (DST), and
 Tile Low-Rank (TLR) computation, and,
3.  Predicting unknown observations using the estimated theta parameters of the Maximum likelihood function (i.e., variance, smoothness, range).
