module switch PrgEnv-cray PrgEnv-gnu
module load cmake
module unload  cray-mpich
module load cray-mpich/7.7.6
export LC_ALL=en_US.UTF-8
export CRAYPE_LINK_TYPE=dynamic
==================================================
export SETUP_DIR=$PWD
rm -rf exageostatr
cd $SETUP_DIR
if [ ! -d "nlopt-2.4.2" ]; then
        wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz
        tar -zxvf nlopt-2.4.2.tar.gz
fi
cd nlopt-2.4.2
[[ -d nlopt_install ]] || mkdir nlopt_install
CC=gcc ./configure --prefix=$PWD/nlopt_install/ --enable-shared --without-guile
make -j
make -j install
export NLOPTROOT=$PWD
export PKG_CONFIG_PATH=$NLOPTROOT/nlopt_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$NLOPTROOT/nlopt_install/lib:$LD_LIBRARY_PATH
#================================
cd $SETUP_DIR
if [ ! -d "gsl-2.4" ]; then
        wget https://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz
        tar -zxvf gsl-2.4.tar.gz
fi
cd gsl-2.4
[[ -d gsl_install ]] || mkdir gsl_install
CC=gcc ./configure --prefix=$PWD/gsl_install/
make -j
make -j install
GSLROOT=$PWD
export PKG_CONFIG_PATH=$GSLROOT/gsl_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$GSLROOT/gsl_install/lib:$LD_LIBRARY_PATH
================================
cd $SETUP_DIR
if [  ! -d "hwloc-1.11.8" ]; then
        wget https://download.open-mpi.org/release/hwloc/v2.0/hwloc-1.11.8.tar.gz
        tar -zxvf hwloc-1.11.8.tar.gz
fi
cd hwloc-1.11.8
[[ -d hwloc_install ]] || mkdir hwloc_install
CC=cc CXX=CC ./configure --prefix=$PWD/hwloc_install --disable-libxml2 -disable-pci --enable-shared=no

make -j
make -j install
export HWLOCROOT=$PWD
export PKG_CONFIG_PATH=$HWLOCROOT/hwloc_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$HWLOCROOT/hwloc_install/lib:$LD_LIBRARY_PATH
================================
cd $SETUP_DIR
if [ ! -d "starpu-1.2.9" ]; then
wget http://starpu.gforge.inria.fr/files/starpu-1.2.9/starpu-1.2.9.tar.gz        
tar -zxvf starpu-1.2.9.tar.gz
fi
cd starpu-1.2.9
[[ -d starpu_install ]] || mkdir starpu_install
CC=cc CXX=CC FC=ftn ./configure --prefix=$PWD/starpu_install/ --disable-cuda --disable-opencl --disable-shared --disable-build-doc --disable-export-dynamic   --with-mpicc=/opt/cray/pe/craype/2.6.3/bin/cc
make -j
make -j  install
export STARPUROOT=$PWD
export PKG_CONFIG_PATH=$STARPUROOT/starpu_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$STARPUROOT/starpu_install/lib:$LD_LIBRARY_PATH
export CPATH=$STARPUROOT/starpu_install/include:$CPATH
#************************************************************************ Install Chameleon - Stars-H - HiCMA
module load intel  #for MKL
cd $SETUP_DIR
# Check if we are already in exageostat repo dir or not.
if git -C $PWD remote -v | grep -q 'https://github.com/ecrc/exageostatr'
then
        # we are, lets go to the top dir (where .git is)
        until test -d $PWD/.git ;
        do
                cd ..
        done;
else
    git clone https://github.com/sc-multivariate-modeling/multivariate-modeling
        cd multivariate-modeling
fi
git pull
git submodule update --init --recursive

export MULTIVARIATEDIR=$PWD/multivariate-modeling
export HICMADIR=$PWD/hicma
export CHAMELEONDIR=$PWD/hicma/chameleon
export STARSHDIR=$PWD/stars-h

## STARS-H
cd $STARSHDIR
rm -rf build
mkdir -p build/install_dir
cd build
cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DMPI=OFF -DOPENMP=OFF -DSTARPU=ON -DBLAS_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}"  -DEXAMPLES=OFF -DTESTING=OFF -DBUILD_SHARED_LIBS=OFF

make -j
make install
export PKG_CONFIG_PATH=$STARSHDIR/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$STARSHDIR/build/install_dir/lib/pkgconfig:$LD_LIBRARY_PATH

## CHAMELEON
cd $CHAMELEONDIR
rm -rf build
mkdir -p build/install_dir
cd build


LDFLAGS=-lrt cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=OFF -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DCHAMELEON_USE_MPI=ON -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_USE_MAGMA=OFF -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON -DCHAMELEON_USE_FXT=OFF -DSTARPU_DIR=$STARPUROOT/starpu_install -DBLAS_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}" -DMORSE_VERBOSE_FIND_PACKAGE=ON -DMPI_C_COMPILER=/opt/cray/pe/craype/2.6.3/bin/cc


make -j # CHAMELEON parallel build seems to be fixed
make install

export PKG_CONFIG_PATH=$CHAMELEONDIR/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$CHAMELEONDIR/build/install_dir/lib/:$LD_LIBRARY_PATH
export CPATH=$CHAMELEONDIR/build/install_dir/include/coreblas:$CPATH
export PATH=$CHAMELEONDIR/build/install_dir/include/coreblas:$PATH
## HICMA
cd $HICMADIR
rm -rf build
mkdir -p build/install_dir
cd build
===============

cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DHICMA_USE_MPI=1 -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=OFF  -DSTARPU_DIR=$STARPUROOT/starpu_install -DBLAS_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core; -lmkl_sequential;-lpthread;-lm;-ldl" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core; -lmkl_sequential;-lpthread; -lm;-ldl" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}" -DMPI_C_COMPILER=/opt/cray/pe/craype/2.6.3/bin/cc

make -j
make install

export PKG_CONFIG_PATH=$HICMADIR/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$HICMADIR/build/install_dir/lib/:$LD_LIBRARY_PATH

$SETUP_DIR

=======================

#Multivariate

export PATH=$CHAMELEONDIR/build/install_dir/include/coreblas:$PATH

cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DEXAGEOSTAT_USE_MPI=1 -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=OFF  -DSTARPU_DIR=$STARPUROOT/starpu_install -DMPI_C_COMPILER=/opt/cray/pe/craype/2.6.3/bin/cc -DEXAGEOSTAT_USE_HICMA=ON -DEXAGEOSTAT_USE_NETCDF=ON -DLAPACKE_INCLUDE_DIRS="-m64;-I${MKLROOT}/include"  -DBLAS_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core; -lmkl_sequential;-lpthread;-lm;-ldl" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core; -lmkl_sequential;-lpthread; -lm;-ldl" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}"
© 2020 GitHub, Inc.
