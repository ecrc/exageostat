export SETUP_DIR=$PWD
module load gcc/4.8.5
module load openblas/0.3.5-nothreads
module load cuda/9.2.148
module load cmake/3.9.2
pause_info(){
    echo "Please press enter key to proceed"; read
    echo "================="
}
    pause_info

echo 'module load gcc/4.8.5' >> $SETUP_DIR/pkg_config.sh
echo 'module load openblas/0.3.5-nothreads' >> $SETUP_DIR/pkg_config.sh
echo 'module load cuda/9.2.148' >> $SETUP_DIR/pkg_config.sh
echo 'module load cmake/3.9.2' >> $SETUP_DIR/pkg_config.sh
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
 ./configure --prefix=$PWD/nlopt_install/ --enable-shared --without-guile
make -j
make -j install
export NLOPTROOT=$PWD
export PKG_CONFIG_PATH=$NLOPTROOT/nlopt_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$NLOPTROOT/nlopt_install/lib:$LD_LIBRARY_PATH
echo 'export PKG_CONFIG_PATH='$NLOPTROOT'/nlopt_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$NLOPTROOT'/nlopt_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
pause_info
======================FXT
========================GSL
module load gsl
========================hwloc
cd $SETUP_DIR
if [  ! -d "hwloc-1.11.5" ]; then
        wget https://download.open-mpi.org/release/hwloc/v1.11/hwloc-1.11.5.tar.gz
        tar -zxvf hwloc-1.11.5.tar.gz
fi
cd hwloc-1.11.5
[[ -d hwloc_install ]] || mkdir hwloc_install

 ./configure --prefix=$PWD/hwloc_install --disable-libxml2 -disable-pci --enable-shared=yes

make -j
make -j install
export HWLOCROOT=$PWD
export PKG_CONFIG_PATH=$HWLOCROOT/hwloc_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$HWLOCROOT/hwloc_install/lib:$LD_LIBRARY_PATH

echo 'export PKG_CONFIG_PATH='$HWLOCROOT'/hwloc_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$HWLOCROOT'/hwloc_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
pause_info

#module load hwloc
==========================starpu
cd $SETUP_DIR
if [ ! -d "starpu-1.2.8" ]; then
        wget http://starpu.gforge.inria.fr/files/starpu-1.2.8/starpu-1.2.8.tar.gz
        tar -zxvf starpu-1.2.8.tar.gz
fi
cd starpu-1.2.8
rm -rf starpu_install
[[ -d starpu_install ]] || mkdir starpu_install

LDFLAGS=-L/sw/summit/cuda/9.2.148/lib64 FC=gfortran CC=gcc CXX=g++ ./configure --prefix=$PWD/starpu_install --disable-opencl --enable-cuda --enable-maxcudadev=6  --with-cuda-dir=/sw/summit/cuda/9.2.148 --with-cuda-include-dir=/sw/summit/cuda/9.2.148/include --with-cuda-lib-dir=/sw/summit/cuda/9.2.148/lib64 --disable-build-doc --disable-export-dynamic --disable-mpi-check  --disable-full-gdb-information

make -j
make -j  install
export STARPUROOT=$PWD
export PKG_CONFIG_PATH=$STARPUROOT/starpu_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$STARPUROOT/starpu_install/lib:$LD_LIBRARY_PATH
export CPATH=$STARPUROOT/starpu_install/include/starpu/1.2/:$CPATH

echo 'export PKG_CONFIG_PATH='$STARPUROOT'/starpu_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$STARPUROOT'/starpu_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export CPATH='$STARPUROOT'/starpu_install/include/starpu/1.2/:$CPATH' >> $SETUP_DIR/pkg_config.sh
    pause_info
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
export CHAMELEONDIR=$EXAGEOSTATDEVDIR/hicma/chameleon
=================
## CHAMELEON
cd $CHAMELEONDIR
rm -rf build
mkdir -p build/install_dir
cd build
 cmake .. -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=OFF -DCHAMELEON_ENABLE_EXAMPLE=OFF -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DCHAMELEON_USE_MPI=ON -DCHAMELEON_USE_CUDA=ON -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON -DSTARPU_DIR=$STARPUROOT/starpu_install
make -j 4 # CHAMELEON parallel build seems to be fixed
make install

export PKG_CONFIG_PATH=$CHAMELEONDIR/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$CHAMELEONDIR/build/install_dir/lib/:$LD_LIBRARY_PATH
export CPATH=$CHAMELEONDIR/build/install_dir/include/coreblas:$CPATH
export CPATH=/ccs/home/abdullsm/codes/exageostatr/src/hicma/chameleon/build/install_dir/include/:$CPATH

echo 'export PKG_CONFIG_PATH='$CHAMELEONDIR'/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$CHAMELEONDIR'/build/install_dir/lib/:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export CPATH='$CHAMELEONDIR'/build/install_dir/include/coreblas:$CPATH' >> $SETUP_DIR/pkg_config.sh

echo 'export CPATH='$CHAMELEONDIR'/build/install_dir/include:$CPATH' >> $SETUP_DIR/pkg_config.sh
    pause_info

==================
#exageostat
cd ~/develop
git clone https://github.com/ecrc/exageostat-dev
cd exageostat-dev
git submodule update --init --recursive

rm -rf build
mkdir -p build/install_dir
cd build
cmake .. -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=OFF -DEXAGEOSTAT_USE_HICMA=OFF -DEXAGEOSTAT_USE_NETCDF=OFF -DEXAGEOSTAT_USE_MPI=ON -DEXAGEOSTAT_USE_CUDA=OFF -DCMAKE_BUILD_TYPE=Release  -DCFLAGS="/sw/summit/cuda/9.2.148/include"


make  && make install