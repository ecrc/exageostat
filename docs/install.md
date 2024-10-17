Installation
============

## Content
1. ExaGeoStat cmake flags
2. Base Installation
3. Dependencies Installations

### 1. ExaGeoStat cmake flags
```-DCMAKE_INSTALL_PREFIX=path```       : Set the installation path for ExaGeoStat. 

```DEXAGEOSTAT_USE_HICMA=ON/OFF```      : Enables ExaGeoStat to use HiCMA.

```DEXAGEOSTAT_USE_CHAMELEON=ON/OFF```  : Enables ExaGeoStat to use Chameleon.

```-DEXAGEOSTAT_SCHED_STARPU=ON/OFF```  : Enables ExaGeoStat to use StarPu.

```DEXAGEOSTAT_USE_NETCDF=ON```         : Enables ExaGeoStat to use NetCDF

```DEXAGEOSTAT_USE_MPI=ON/OFF```        : Enables ExaGeoStat to use MPI.

```EXAGEOSTAT_USE_CUDA=ON/OFF```        : Enables ExaGeoStat to use Cuda.

```DMPI_VALUE=ON/OFF```                 : Set MPI value for other dependencies.

```DCUDA_VALUE=ON/OFF```                : Set cuda value for other dependencies.

```-DEXAGEOSTAT_PACKAGE=ON/OFF```       : Enable a packaging system for distribution

```DEXAGEOSTAT_INSTALL_DEPS=ON/OFF```   : Set the automatic installations of ExaGeoStat's dependencies

```DBUILD_SHARED_LIBS=ON/OFF```         : Enables building the shared libs.

```DEXAGEOSTAT_EXAMPLES=ON/OFF```       : Enables ExaGeoStat's examples.

### 2. Base Installation
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

        mkdir -p build/installdir && cd build

5.  Use CMake to get all the dependencies

        cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir  -DEXAGEOSTAT_SCHED_STARPU=ON   -DEXAGEOSTAT_USE_NETCDF=ON -DEXAGEOSTAT_USE_HICMA=ON -DEXAGEOSTAT_USE_CHAMELEON=ON -DEXAGEOSTAT_INSTALL_DEPS=ON -DBUILD_SHARED_LIBS=ON

6.  Build EXAGEOSTAT

        make -j

7.  Build local documentation (optional)

        make docs

8.  Install EXAGEOSTAT

        make install

9. Add line

        export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

   to your .bashrc file.

Now you can use pkg-config executable to collect compiler and linker flags for
EXAGEOSTAT.

### 3. Dependencies Installation

#### NetCDF
```commandline
    echo "Installing hdf5 library for netcdf4..."
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.gz -O - | tar -zx
    echo "Building hdf5..."
    cd hdf5*
    ./configure --prefix="$SETUP_DIR" --with-zlib="$SETUP_DIR" --enable-hl
    make -j || { echo 'hdf5 installation failed' ; exit 1; }
    make -j install
    echo "hdf5 Done..."
    
    export HDF5=$SETUP_DIR
    export LD_LIBRARY_PATH=$SETUP_DIR/lib:$LD_LIBRARY_PATH
    
    echo "Installing NetCDF..."
    cd "$SETUP_DIR" || exit
    wget https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.7.4.tar.gz -O - | tar -zx
    echo "Building NetCDF..."
    cd netcdf* || exit
    
    export CPPFLAGS=-I$SETUP_DIR/include
    export LDFLAGS=-L$SETUP_DIR/lib
    
    ./configure --disable-dap --prefix="$SETUP_DIR"
    make -j || { echo 'Netcdf installation failed' ; exit 1; }
    make -j install
    
    export PKG_CONFIG_PATH=$SETUP_DIR/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$SETUP_DIR/lib:$LD_LIBRARY_PATH
```
#### NLOPT
```commandline
    echo "Installing NLOPT..."
    wget https://github.com/stevengj/nlopt/archive/v2.7.0.tar.gz -O - | tar -zx
    echo "Building NLOPT..."
    cd nlopt* || exit
    mkdir -p build
    cd build || exit
    rm -rf ./CMake*
    cmake  -DCMAKE_C_FLAGS=-fcommon -DCMAKE_C_FLAGS=-fPIC -DCMAKE_BUILD_TYPE="Release" \
    -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w" -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX="$SETUP_DIR"  ..
    make -j || { echo 'NLOPT installation failed' ; exit 1; }
    make -j install
    
    export PKG_CONFIG_PATH=$SETUP_DIR/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$SETUP_DIR/lib:$LD_LIBRARY_PATH
    export NLOPT_LIBRARY_DIRS=$SETUP_DIR/lib:$NLOPT_LIBRARY_DIRS
    export NLOPT_INCLUDE_DIRS=$SETUP_DIR/include:NLOPT_INCLUDE_DIRS
```
#### HWLOC
```commandline
    echo "Installing HWLOC..."
    wget https://download.open-mpi.org/release/hwloc/v2.4/hwloc-2.4.0.tar.gz -O - | tar -zx
    echo "Building HWLOC..."
    cd hwloc* || exit
    ./configure --prefix="$SETUP_DIR"
    make -j || { echo 'HWLOC installation failed' ; exit 1; }
    make -j install
    
    export PKG_CONFIG_PATH=$SETUP_DIR/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$SETUP_DIR/lib:$LD_LIBRARY_PATH
```

#### GSL
```commandline
    echo "Installing GSL..."
    wget https://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz -O - | tar -zx
    echo "Building GSL..."
    cd gsl* || exit
    ./configure --prefix="$SETUP_DIR"
    make -j || { echo 'GSL installation failed' ; exit 1; }
    make -j install
    
    export PKG_CONFIG_PATH=$SETUP_DIR/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$SETUP_DIR/lib:$LD_LIBRARY_PATH
```

#### STARPU
```commandline
    echo "Installing StarPu..."
    cd "$SETUP_DIR" || exit
    wget https://files.inria.fr/starpu/starpu-1.3.9/starpu-1.3.9.tar.gz -O - | tar -zx
    echo "Building StarPu... with CUDA $CUDA_VALUE AND MPI $MPI_VALUE"
    cd starpu* || exit
    ./configure --disable-starpufft --enable-cuda  --disable-mpi --disable-opencl --prefix="$SETUP_DIR"   --disable-starpu-top --disable-starpufft --disable-build-doc --disable-starpufft-examples   --disable-fortran --with-perf-model-dir="$SETUP_DIR" --disable-fstack-protector-all --disable-gcc-extensions
    make -j || { echo 'STARPU installation failed' ; exit 1; }
    make -j install
    export PKG_CONFIG_PATH=$SETUP_DIR/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$SETUP_DIR/lib:$LD_LIBRARY_PATH
```
