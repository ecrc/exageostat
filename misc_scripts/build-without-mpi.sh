#!/bin/bash -le
module load ecrc-extras
module load mkl/2018-initial
module load gcc/5.5.0
module load cmake/3.9.6
module load hwloc/1.11.8-gcc-5.5.0
#module load openmpi/3.0.0-gcc-5.5.0
module load fxt/0.3.7-gcc-5.5.0
module load starpu/1.2.6-gcc-5.5.0-mkl-openmpi-3.0.0-fxt
module load gsl/2.4-gcc-5.5.0
module load nlopt/2.4.2-gcc-5.5.0
module load r/3.2.3
module load hdf5/1.10.1-gcc-5.5.0
module load netcdf/4.5.0-gcc-5.5.0


git config --global credential.helper 'cache --timeout=36000'


# BASH verbose mode
#set -x
SETUP_DIR="/home/abdullsm/codes"

echo 'The installation directory is '$SETUP_DIR
echo 'The mkl root directory is '$MKL_DIR

#***************************************************************************** clean bashrc from previous installation
sed -i '/## EXAGEOSTAT-INSTALLATION-BEGIN/,/## EXAGEOSTAT-INSTALLATION-END/d'  ~/.bashrc
#************************************************************************ Install Chameleon - Stars-H - HiCMA 
cd $SETUP_DIR
# Check if we are already in exageostat repo dir or not.
[[ -d trash ]] || mkdir trash
DATE=`date +%Y-%m-%d`
if [ -d exageostat-dev ]; then
	rm -rf ./trash/$DATE/
	mv ./exageostat-dev ./trash/$DATE/
fi
if [ -d hicma ]; then
	rm -rf ./trash/$DATE/
        mv ./hicma ./trash/$DATE/
fi
if [ -d chameleon ]; then
	rm -rf ./trash/$DATE/
        mv ./chameleon ./trash/$DATE/
fi
if [ -d stars-h ]; then
	rm -rf ./trash/$DATE/
        mv ./stars-h ./trash/$DATE/
fi

git clone https://github.com/ecrc/exageostat-dev
cd exageostat-dev
git submodule update --init --recursive
EXAGEOSTATROOT=$PWD
#Copy hicma, and stars-H to set_up directory
mv stars-h $SETUP_DIR
mv hicma  $SETUP_DIR
############################# Chameleon Installation
cd $SETUP_DIR
cd hicma
cd chameleon
CHAMELEONROOT=$PWD
echo $CHAMELEONROOT
mkdir -p build/installdir
cd build
CC=gcc cmake .. -DCMAKE_BUILD_TYPE=Debug -DCHAMELEON_USE_MPI=OFF -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DBUILD_SHARED_LIBS=ON
make -j
make install
export PKG_CONFIG_PATH=$CHAMELEONROOT/build/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$CHAMELEONROOT/build/installdir/lib:$LD_LIBRARY_PATH
############################# Stars-H Installation
cd $SETUP_DIR
cd stars-h
STARSHROOT=$PWD
mkdir -p build/installdir
cd build
CC=gcc cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir/ -DCMAKE_C_FLAGS=-fPIC  -DMPI=OFF -DOPENMP=OFF -DSTARPU=OFF
make -j
make install
export PKG_CONFIG_PATH=$STARSHROOT/build/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$STARSHROOT/build/installdir/lib:$LD_LIBRARY_PATH
############################# HiCMA Installation
cd $SETUP_DIR
cd hicma
HICMAROOT=$PWD
mkdir -p build/installdir
cd build
CC=gcc cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DHICMA_USE_MPI=OFF -DCMAKE_C_FLAGS=-fPIC
make -j
make install
export PKG_CONFIG_PATH=$HICMAROOT/build/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$HICMAROOT/build/installdir/lib:$LD_LIBRARY_PATH
#***************************************************************************** edit bashrc file
echo '## EXAGEOSTAT-INSTALLATION-BEGIN' >> ~/.bashrc
echo 'module load ecrc-extras' >> ~/.bashrc
echo 'module load mkl/2018-initial' >> ~/.bashrc
echo 'module load gcc/5.5.0' >> ~/.bashrc
echo 'module load cmake/3.9.6' >> ~/.bashrc
echo 'module load hwloc/1.11.8-gcc-5.5.0' >> ~/.bashrc
#echo 'module load openmpi/3.0.0-gcc-5.5.0' >> ~/.bashrc
echo 'module load fxt/0.3.7-gcc-5.5.0' >> ~/.bashrc
echo 'module load starpu/1.2.6-gcc-5.5.0-mkl-openmpi-3.0.0-fxt' >> ~/.bashrc
echo 'module load gsl/2.4-gcc-5.5.0' >> ~/.bashrc
echo 'module load nlopt/2.4.2-gcc-5.5.0' >> ~/.bashrc
echo 'module load r/3.2.3' >> ~/.bashrc
echo 'module load hdf5/1.10.1-gcc-5.5.0' >> ~/.bashrc
echo 'module load netcdf/4.5.0-gcc-5.5.0' >> ~/.bashrc

#MKL
#echo '. '$MKLROOT'/bin/mklvars.sh intel64' >> ~/.bashrc
#echo 'export MKLROOT='$MKL_DIR >> ~/.bashrc
#echo 'export LD_PRELOAD='$MKL_DIR'/lib/intel64/libmkl_core.so:'$MKL_DIR'/lib/intel64/libmkl_sequential.so' >> ~/.bashrc

#starpu
echo 'export PKG_CONFIG_PATH='$STARPUROOT'/starpu_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH='$STARPUROOT'/starpu_install/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
#CHAMELEON
echo 'export PKG_CONFIG_PATH='$CHAMELEONROOT'/build/installdir/lib/pkgconfig:$PKG_CONFIG_PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH='$CHAMELEONROOT'/build/installdir/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
#HICMA
echo 'export PKG_CONFIG_PATH='$HICMAROOT'/build/installdir/lib/pkgconfig:$PKG_CONFIG_PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH='$HICMAROOT'/build/installdir/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
echo 'export CPATH='$HICMAROOT'/build/installdir/include:$CPATH' >> ~/.bashrc
#STARS-H
echo 'export PKG_CONFIG_PATH='$STARSHROOT'/build/installdir/lib/pkgconfig:$PKG_CONFIG_PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH='$STARSHROOT'/build/installdir/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
echo 'export CPATH='$STARSHROOT'/build/installdir/include:$CPATH' >> ~/.bashrc
#end
echo '## EXAGEOSTAT-INSTALLATION-END' >> ~/.bashrc
##################################################################################
