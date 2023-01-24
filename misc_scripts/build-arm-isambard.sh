export SETUP_DIR=$PWD
module load cdt/19.08
module swap PrgEnv-cray/6.0.5 PrgEnv-allinea

ARMPLLIBS="-larmpl"
ARMPLROOT="/opt/allinea/19.2.0.0/opt/arm/armpl-19.2.0_ThunderX2CN99_SUSE-12_arm-hpc-compiler_19.2_aarch64-linux" >> $SETUP_DIR/pkg_config.sh
export LD_LIBRARY_PATH="/opt/allinea/19.2.0.0/opt/arm/armpl-19.2.0_ThunderX2CN99_SUSE-12_arm-hpc-compiler_19.2_aarch64-linux/lib":$LD_LIBRARY_PATH


echo 'module load cdt/19.08' >> $SETUP_DIR/pkg_config.sh
echo 'module swap PrgEnv-cray PrgEnv-allinea' >> $SETUP_DIR/pkg_config.sh
echo 'ARMPLLIBS="-larmpl"' >> $SETUP_DIR/pkg_config.sh
echo 'ARMPLROOT="/opt/allinea/19.2.0.0/opt/arm/armpl-19.2.0_ThunderX2CN99_SUSE-12_arm-hpc-compiler_19.2_aarch64-linux"' >> $SETUP_DIR/pkg_config.sh


================== NLOPT
export SETUP_DIR=$PWD
rm -rf exageostatr
cd $SETUP_DIR
if [ ! -d "nlopt-2.4.2" ]; then
        wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz
        tar -zxvf nlopt-2.4.2.tar.gz
fi
cd nlopt-2.4.2
[[ -d nlopt_install ]] || mkdir nlopt_install
CC=cc CXX=CC ./configure --prefix=$PWD/nlopt_install/ --enable-shared --without-guile
make -j
make -j install
export NLOPTROOT=$PWD
export PKG_CONFIG_PATH=$NLOPTROOT/nlopt_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$NLOPTROOT/nlopt_install/lib:$LD_LIBRARY_PATH
echo 'export PKG_CONFIG_PATH='$NLOPTROOT'/nlopt_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$NLOPTROOT'/nlopt_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
======================FXT
========================GSL
cd $SETUP_DIR
if [ ! -d "gsl-2.4" ]; then
        wget https://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz
        tar -zxvf gsl-2.4.tar.gz
fi
cd gsl-2.4
[[ -d gsl_install ]] || mkdir gsl_install
CC=cc CXX=CC ./configure --prefix=$PWD/gsl_install/
make -j
make -j install
export GSLROOT=$PWD
export PKG_CONFIG_PATH=$GSLROOT/gsl_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$GSLROOT/gsl_install/lib:$LD_LIBRARY_PATH
echo 'export PKG_CONFIG_PATH='$GSLROOT'/gsl_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$GSLROOT'/gsl_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
========================hwloc
cd $SETUP_DIR
if [  ! -d "hwloc-1.11.5" ]; then
        wget https://download.open-mpi.org/release/hwloc/v2.0/hwloc-1.11.5.tar.gz
        tar -zxvf hwloc-1.11.5.tar.gz
fi
cd hwloc-1.11.5
[[ -d hwloc_install ]] || mkdir hwloc_install
CC=cc CXX=CC ./configure --prefix=$PWD/hwloc_install --disable-libxml2 -disable-pci --enable-shared=yes

make -j
make -j install
export HWLOCROOT=$PWD
export PKG_CONFIG_PATH=$HWLOCROOT/hwloc_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$HWLOCROOT/hwloc_install/lib:$LD_LIBRARY_PATH

echo 'export PKG_CONFIG_PATH='$HWLOCROOT'/hwloc_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$HWLOCROOT'/hwloc_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
==========================starpu
cd $SETUP_DIR
if [ ! -d "starpu-1.2.6" ]; then
        wget http://starpu.gforge.inria.fr/files/starpu-1.2.6/starpu-1.2.6.tar.gz
        tar -zxvf starpu-1.2.6.tar.gz
fi
cd starpu-1.2.6
[[ -d starpu_install ]] || mkdir starpu_install
CC=cc CXX=CC ./configure --prefix=$PWD/starpu_install/ --disable-cuda --disable-opencl  --disable-build-doc --disable-export-dynamic --disable-mpi-check --with-mpicc=`which cc`
make -j
make -j  install
export STARPUROOT=$PWD
export PKG_CONFIG_PATH=$STARPUROOT/starpu_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$STARPUROOT/starpu_install/lib:$LD_LIBRARY_PATH
export CPATH=$STARPUROOT/starpu_install/include/starpu/1.2/:$CPATH

echo 'export PKG_CONFIG_PATH='$STARPUROOT'/starpu_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$STARPUROOT'/starpu_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export CPATH='$STARPUROOT'/starpu_install/include/starpu/1.2/:$CPATH' >> $SETUP_DIR/pkg_config.sh



=================================
cd $SETUP_DIR
rm -rf exageostatr
# Check if we are already in exageostat repo dir or not.
if git -C $PWD remote -v | grep -q 'https://github.com/ecrc/exageostatr'
then
        # we are, lets go to the top dir (where .git is)
        until test -d $PWD/.git ;
        do
                cd ..
        done;
else
        git clone https://github.com/ecrc/exageostatr
        cd exageostatr
fi
git pull
git submodule update --init --recursive

export EXAGEOSTATDEVDIR=$PWD/src
cd $EXAGEOSTATDEVDIR
export HICMADIR=$EXAGEOSTATDEVDIR/hicma
export CHAMELEONDIR=$EXAGEOSTATDEVDIR/hicma/chameleon
export STARSHDIR=$EXAGEOSTATDEVDIR/stars-h


==========================================
## CHAMELEON
cd $CHAMELEONDIR
rm -rf build
mkdir -p build/install_dir
cd build


cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=OFF -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=OFF -DCHAMELEON_ENABLE_TIMING=ON -DCHAMELEON_USE_MPI=ON -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON -DSTARPU_DIR=$STARPUROOT/starpu_install -DBLAS_LIBRARIES="-Wl,--no-as-needed;-L${ARMPLROOT}/lib;${ARMPLLIBS};-lpthread;-lm;-ldl" -DBLAS_COMPILER_FLAGS="-I${ARMPLROOT}/include" -DLAPACK_LIBRARIES="-Wl,--no-as-needed;-L${ARMPLROOT}/lib;${ARMPLLIBS} -lpthread;-lm;-ldl" -DCBLAS_DIR="${ARMPLROOT}" -DLAPACKE_DIR="${ARMPLROOT}"   -DMPI_C_COMPILER=`which cc`


make -j 4 # CHAMELEON parallel build seems to be fixed
make install

export PKG_CONFIG_PATH=$CHAMELEONDIR/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH

echo 'export PKG_CONFIG_PATH='$CHAMELEONDIR'/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh

==================
#exageostat
cd ~/develop
git clone https://github.com/ecrc/exageostat-dev
cd exageostat-dev
rm -rf build
mkdir -p build/install_dir
cd build
module load cray-netcdf-hdf5parallel/4.6.3.0

cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=OFF -DSTARPU_DIR=$STARPUROOT/starpu_install -DBLAS_LIBRARIES="-Wl,--no-as-needed;-L${ARMPLROOT}/lib;${ARMPLLIBS};-lpthread;-lm;-ldl" -DBLAS_COMPILER_FLAGS="-I${ARMPLROOT}/include" -DLAPACK_LIBRARIES="-Wl,--no-as-needed;-L${ARMPLROOT}/lib;${ARMPLLIBS} -lpthread;-lm;-ldl" -DCBLAS_DIR="${ARMPLROOT}" -DLAPACKE_DIR="${ARMPLROOT}"   -DMPI_C_COMPILER=`which cc`  -DEXAGEOSTAT_USE_HICMA=OFF -DEXAGEOSTAT_USE_NETCDF=OFF -DEXAGEOSTAT_USE_MPI=1

make -j 4 && make install