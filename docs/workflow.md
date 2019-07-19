Workflow of ExaGeoStat
======================
The Exascale Geo-Statistics project (ExaGeoStat) is a parallel high performance
unified framework for computational geostatistics on manycore systems. The
software project aims at optimizing the likelihood function for a given spatial data
in order to provide an efficient way for predicting missing observations in the context of
climate/weather forecasting applications. The framework has three main components:


1. Generate synthetic datasets in the case of testing mode. The synthetic dataset
   includes 2D spatial locations and measurements vector. The user can enable testing
   mode by using (--test) input. Furthermore, real datasets can be used with
   ExaGeoStat, when passing input files names for both spatial locations and measurements.

2. An iterative procedure is performed over the given data (i.e., synthetic or real dataset)
   to estimate the parameter vector the maximize the likelihood function for the given
   data. The average number of iterations depends on the lower and upper bounds of the
   optimization algorithm given by the user. Three computation options are available here,
   Exact, Diagonal Super Tile (DST), and Tile Low-Rank (TLR)

3. The estimated parameter vector (i.e., variance, smoothness, and range) of the Matern matrix kernel can be used
   to predict measurements on missing spatial locations. The ExaGeoStat's prediction stage
   is able to predict the measurements from the MLE parameter vectors. This stage is also running 
   in two different modes: test and real. For testing, unobserved spatial locations are randomly chosen from
   the original spatial locations and their measurements are marked as unknown. In this
   case, ExaGeoStat is able to compute the Mean Square Error (MSE) by comparing observed
   measurements with predicted measurements. For real cases, the user provides the program
   with the new spatial locations and the program predict the measurements for those new locations.
