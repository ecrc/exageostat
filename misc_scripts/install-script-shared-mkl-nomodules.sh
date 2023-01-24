#!/bin/bash -x


# variables
export BASEDIR=$PWD
export SETUP_DIR="/home/qadirga/codes"
export TMPDIR=/tmp/_$$

mkdir -p $TMPDIR

#SETUP_DIR=$1
if [ -z "$SETUP_DIR" ]; then
# Use RLIBS for setup dir
arr=(`Rscript -e '.libPaths()' | gawk '{printf "%s ",$2}'`)
for i in ${!arr[*]};
do
    dir=`echo ${arr[$i]}|tr -d \"`
    if [ -d "$dir" ] && [ -w "$dir" ]
    then
        SETUP_DIR="$dir/exageostat"
        break
    fi
done
fi
mkdir -p $SETUP_DIR

if [ -z "$SETUP_DIR" ]
then
    echo "Check your .libPaths() in R. Could not find a writable directory."
    exit 1;
fi

if [ -n "$MKLROOT" ] && [ -d "$MKLROOT" ]; then
    echo "mkl_dir directory exists!"
    echo "Great... continue set-up"
else
    echo "MKLROOT Directory does not exist!... Please define and export MKLROOT variable"
    exit 1
fi
PREFIX=$SETUP_DIR


echo 'The installation directory is '$SETUP_DIR
echo 'The mkl root directory is '$MKLROOT

############################## Check OS
echo "Finding the current os type"
echo
osType=$(uname)
case "$osType" in 
	"Darwin")
	{
	  echo "Running on Mac OSX."
	  CURRENT_OS="OSX"
	} ;;    
	"Linux")
	{
	  echo "Running on LINUX."
	  CURRENT_OS="LINUX"
	} ;;
       *) 
       {
         echo "Unsupported OS, exiting"
         exit
       } ;;

esac

#################################################
export MKLROOT=$MKLROOT
. $MKLROOT/bin/mklvars.sh intel64



export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH
rpaths="-Wl,-rpath=$PREFIX/lib -Wl,-rpath=$PREFIX/libs -L$PREFIX/lib "
echo "LDFLAGS += $rpaths " >> src/Makefile

#*****************************************************************************
set -e

if [ $CURRENT_OS == "LINUX" ]
then
    export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
else
    export DYLD_LIBRARY_PATH=$PREFIX/lib:$DYLD_LIBRARY_PATH
fi

#*****************************************************************************install Nlopt
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

#========================hwloc
cd $SETUP_DIR
if [  ! -d "hwloc-2.0.2" ]; then
        wget https://download.open-mpi.org/release/hwloc/v2.0/hwloc-2.0.2.tar.gz
        tar -zxvf hwloc-2.0.2.tar.gz
fi
cd hwloc-2.0.2
[[ -d hwloc_install ]] || mkdir hwloc_install
./configure --prefix=$PWD/hwloc_install --disable-libxml2 -disable-pci --enable-shared=yes

make -j
make -j install
export HWLOCROOT=$PWD
export PKG_CONFIG_PATH=$HWLOCROOT/hwloc_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$HWLOCROOT/hwloc_install/lib:$LD_LIBRARY_PATH
echo 'export PKG_CONFIG_PATH='$HWLOCROOT'/hwloc_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$HWLOCROOT'/hwloc_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
#================================
cd $SETUP_DIR
if [ ! -d "gsl-2.4" ]; then
        wget https://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz
        tar -zxvf gsl-2.4.tar.gz
fi
cd gsl-2.4
[[ -d gsl_install ]] || mkdir gsl_install
 ./configure --prefix=$PWD/gsl_install/
make -j
make -j install
export GSLROOT=$PWD
export PKG_CONFIG_PATH=$GSLROOT/gsl_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$GSLROOT/gsl_install/lib:$LD_LIBRARY_PATH

echo 'export PKG_CONFIG_PATH='$GSLROOT'/gsl_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$GSLROOT'/gsl_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
#================================
#================================
cd $SETUP_DIR
if [ ! -d "starpu-1.2.5" ]; then
        wget http://starpu.gforge.inria.fr/files/starpu-1.2.5/starpu-1.2.5.tar.gz
        tar -zxvf starpu-1.2.5.tar.gz
fi
cd starpu-1.2.5
[[ -d starpu_install ]] || mkdir starpu_install
./configure --prefix=$SETUP_DIR/starpu-1.2.5/starpu_install  -disable-cuda --disable-opencl --enable-shared --disable-build-doc --disable-export-dynamic --disable-mpi-check
make -j
make -j  install
export STARPUROOT=$PWD
export PKG_CONFIG_PATH=$STARPUROOT/starpu_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$STARPUROOT/starpu_install/lib:$LD_LIBRARY_PATH
export CPATH=$STARPUROOT/starpu_install/include/starpu/1.2:$CPATH
echo 'export PKG_CONFIG_PATH='$STARPUROOT'/starpu_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$STARPUROOT'/starpu_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export CPATH='$STARPUROOT'/starpu_install/include/starpu/1.2:$CPATH' >> $SETUP_DIR/pkg_config.sh

#************************************************************************ Install Chameleon - Stars-H - HiCMA
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
        git clone https://github.com/ecrc/exageostatr
        cd exageostatr
fi
git pull
git submodule update --init --recursive

export EXAGEOSTATDEVDIR=$PWD/src
cd $EXAGEOSTATDEVDIR
rm -rf hicma
rm -rf stars-h
git clone https://github.com/SAbdulah/hicma
cd hicma
git submodule update --init --recursive
cd ..
git clone https://github.com/SAbdulah/stars-h
cd stars-h
git submodule update --init --recursive
cd ..
export HICMADIR=$EXAGEOSTATDEVDIR/hicma
export CHAMELEONDIR=$EXAGEOSTATDEVDIR/hicma/chameleon
export STARSHDIR=$EXAGEOSTATDEVDIR/stars-h

## STARS-H
cd $STARSHDIR
rm -rf build
mkdir -p build/install_dir


cd build

cmake ..  -DCMAKE_INSTALL_PREFIX=$STARSHDIR/build/install_dir -DMPI=OFF -DOPENMP=OFF -DSTARPU=OFF -DCMAKE_C_FLAGS="-fPIC"

make -j
make install

export PKG_CONFIG_PATH=$STARSHDIR/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$STARSHDIR/build/install_dir/lib:$LD_LIBRARY_PATH

echo 'export PKG_CONFIG_PATH='$STARSHDIR'/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export PKG_CONFIG_PATH='$STARSHDIR'/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH' >>  $SETUP_DIR/pkg_config.sh

## CHAMELEON
cd $CHAMELEONDIR
rm -rf build
mkdir -p build/install_dir
cd build


 cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/install_dir  -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DCHAMELEON_USE_MPI=OFF -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_USE_MAGMA=OFF -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON -DCHAMELEON_USE_FXT=OFF -DSTARPU_DIR=$STARPUROOT/starpu_install  -DCHAMELEON_VERBOSE_FIND_PACKAGE=ON

make -j # CHAMELEON parallel build seems to be fixed
make install

export PKG_CONFIG_PATH=$CHAMELEONDIR/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$CHAMELEONDIR/build/install_dir/lib:$LD_LIBRARY_PATH
export CPATH=$CHAMELEONDIR/build/install_dir/include/coreblas:$CPATH


echo 'export PKG_CONFIG_PATH='$CHAMELEONDIR'/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH'  >>  $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$CHAMELEONDIR'/build/install_dir/lib:$LD_LIBRARY_PATH'  >>  $SETUP_DIR/pkg_config.sh
echo 'export CPATH='$CHAMELEONDIR'/build/install_dir/include/coreblas:$CPATH'  >>  $SETUP_DIR/pkg_config.sh

#HICMA
cd $HICMADIR
rm -rf build
mkdir -p build/install_dir
cd build
#===============
cmake ..  -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON -DHICMA_USE_MPI=OFF 
 
make -j
make install

export PKG_CONFIG_PATH=$HICMADIR/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$HICMADIR/build/install_dir/lib:$LD_LIBRARY_PATH


echo 'export PKG_CONFIG_PATH='$HICMADIR'/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH' >>  $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$HICMADIR'/build/install_dir/lib:$LD_LIBRARY_PATH' >>  $SETUP_DIR/pkg_config.sh

