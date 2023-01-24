module load gcc/6.4.0
module load cuda/10.1.105
module load openmpi/4.0.1/gnu6.4.0_cuda10.1.105
git clone https://github.com/ecrc/exageostat
export EXAGEOSTATDEVDIR=$PWD/exageostat
mkdir codes
SETUP_DIR=$(pwd)/codes
CUDAVALUE="ON"
MPIVALUE="ON"
BUILD_DEPENDENCIES="true"

cd exageostat
git pull
git submodule update --init --recursive

cd $SETUP_DIR

mv ../exageostat/hicma ./
mv ../exageostat/stars-h ./
export HICMADIR=$PWD/hicma
export CHAMELEONDIR=$PWD/hicma/chameleon
export STARSHDIR=$PWD/stars-h

# openblas
if pkg-config --exists --atleast-version=2 openblas
then
    _LOCATION=`pkg-config --variable=prefix openblas`
    echo "gsl FOUND in [$_LOCATION]"
else
    if [ "$BUILD_DEPENDENCIES" == "true" ]
    then
        echo "Building OpenBLAS..."
        cd $SETUP_DIR
        wget https://github.com/xianyi/OpenBLAS/archive/v0.3.7.tar.gz -O - | tar -zx
        cd OpenBLAS-0.3.7
    	mkdir openblas_install
        make PREFIX=$PWD/openblas_install >/dev/null|| make PREFIX=$PWD/openblas_install || { echo 'OpenBLAS installation failed' ; exit 1; }
        make install PREFIX=$PWD/openblas_install
        export OPENBLASROOT=$PWD
        export PKG_CONFIG_PATH=$OPENBLASROOT/openblas_install/lib/pkgconfig:$PKG_CONFIG_PATH
        export LD_LIBRARY_PATH=$OPENBLASROOT/openblas_install/lib:$LD_LIBRARY_PATH
        echo 'export PKG_CONFIG_PATH='$OPENBLASROOT'/openblas_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
        echo 'export LD_LIBRARY_PATH='$OPENBLASROOT'/openblas_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
    else
        echo "####################"
        echo "openblas NOT FOUND"
        echo "Please download it from: https://github.com/xianyi/OpenBLAS/archive/v0.3.7.tar.gz"
        echo "After installing it, set the proper PKG_CONFIG_PATH variable"
        echo ""
        err=1
    fi
fi

# gsl
if pkg-config --exists --atleast-version=2 gsl
then
    _LOCATION=`pkg-config --variable=prefix gsl`
    echo "gsl FOUND in [$_LOCATION]"
else
    if [ "$BUILD_DEPENDENCIES" == "true" ]
    then
        echo "Building GSL..."
        cd $SETUP_DIR
        wget https://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz -O - | tar -zx
        cd gsl-2.4
    mkdir gsl_install
    ./configure --prefix=$PWD/gsl_install/
    make || make VERBOSE=1 || { echo 'GSL installation failed' ; exit 1; }
    make install
    export GSLROOT=$PWD
    export PKG_CONFIG_PATH=$GSLROOT/gsl_install/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$GSLROOT/gsl_install/lib:$LD_LIBRARY_PATH
    echo 'export PKG_CONFIG_PATH='$GSLROOT'/gsl_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
    echo 'export LD_LIBRARY_PATH='$GSLROOT'/gsl_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
    else
        echo "####################"
        echo "gsl NOT FOUND"
        echo "Please download it from: https://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz"
        echo "After installing it, set the proper PKG_CONFIG_PATH variable"
        echo ""
        err=1
    fi
fi

# nlopt
if pkg-config --exists --atleast-version=2.4 nlopt
then
    _LOCATION=`pkg-config --variable=prefix nlopt`
    echo "nlopt FOUND in [$_LOCATION]"
else
    if [ "$BUILD_DEPENDENCIES" == "true" ]
    then
        echo "Building NLOPT..."
        cd $SETUP_DIR
        wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz -O - | tar -zx
        cd nlopt-2.4.2
    [[ -d nlopt_install ]] || mkdir nlopt_install
        ./configure --enable-shared --without-guile --prefix=$PWD/nlopt_install/
        make || make VERBOSE=1 || { echo 'NLOPT installation failed' ; exit 1; }
        make install
    NLOPTROOT=$PWD
    export PKG_CONFIG_PATH=$NLOPTROOT/nlopt_install/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$NLOPTROOT/nlopt_install/lib:$LD_LIBRARY_PATH
    echo 'export PKG_CONFIG_PATH='$NLOPTROOT'/nlopt_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
    echo 'export LD_LIBRARY_PATH='$NLOPTROOT'/nlopt_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
    else
        echo "####################"
        echo "nlopt NOT FOUND"
        echo "Please download it from: http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz"
        echo "After installing it, set the proper PKG_CONFIG_PATH variable"
        echo ""
        err=1
    fi
fi

# hwloc
if pkg-config --exists hwloc
then
    _LOCATION=`pkg-config --variable=prefix hwloc`
    echo "hwloc FOUND in [$_LOCATION]"
else
    if [ "$BUILD_DEPENDENCIES" == "true" ]
    then
        echo "Building Hwloc..."
        cd $SETUP_DIR
        wget https://download.open-mpi.org/release/hwloc/v2.0/hwloc-2.0.2.tar.gz -O - | tar -zx
	cd hwloc-2.0.2
	[[ -d hwloc_install ]] || mkdir hwloc_install
	CC=gcc ./configure --prefix=$PWD/hwloc_install --disable-libxml2 -disable-pci --enable-shared=yes
       
	 make || make VERBOSE=1 || { echo 'HWLOC installation failed' ; exit 1; }
        make install
    export HWLOCROOT=$PWD
    export PKG_CONFIG_PATH=$HWLOCROOT/hwloc_install/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$HWLOCROOT/hwloc_install/lib:$LD_LIBRARY_PATH
    echo 'export PKG_CONFIG_PATH='$HWLOCROOT'/hwloc_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
    echo 'export LD_LIBRARY_PATH='$HWLOCROOT'/hwloc_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
    else
        echo "####################"
        echo "hwloc NOT FOUND"
        echo "Please download it from: https://www.open-mpi.org/software/hwloc/v1.11/downloads/hwloc-1.11.5.tar.gz"
        echo "After installing it, set the proper PKG_CONFIG_PATH variable"
        echo ""
        err=1
    fi
fi

# StarPU
if pkg-config --exists  starpu-1.2

then
    _LOCATION=`pkg-config --variable=prefix libstarpu`
    echo "StarPU FOUND in [$_LOCATION]"
else
    if [ "$BUILD_DEPENDENCIES" == "true" ]
    then
        echo "Building StarPU..."
    cd $SETUP_DIR
        wget http://starpu.gforge.inria.fr/files/starpu-1.2.8/starpu-1.2.8.tar.gz -O - | tar -zx
        cd starpu-1.2.8
    [[ -d starpu_install ]] || mkdir starpu_install
sed -i 's/MPI_Address/MPI_Get_address/g' mpi/examples/user_datatype/my_interface.c
        if [ "$CUDAVALUE" == "ON" ]; then
                if [ "$MPIVALUE" == "ON" ]; then
                        CFLAGS=-w ./configure --disable-starpufft --enable-cuda --disable-opencl --prefix=$PWD/starpu_install  --enable-blas-lib=none --disable-starpu-top --disable-starpufft --disable-build-doc --disable-starpufft-examples   --disable-fortran --with-perf-model-dir=$PWD/starpu_install 
                else
                        CFLAGS=-w ./configure --disable-starpufft --enable-cuda --disable-opencl --prefix=$PWD/starpu_install  --enable-blas-lib=none --disable-starpu-top --disable-starpufft --disable-build-doc --disable-starpufft-examples --disable-mpi  --disable-fortran --with-perf-model-dir=$PWD/starpu_install 
                fi
        else
                if [ "$MPIVALUE" == "ON" ]; then
                        CFLAGS=-w ./configure --disable-starpufft --disable-cuda --disable-opencl --prefix=$PWD/starpu_install --enable-blas-lib=none --disable-starpu-top --disable-starpufft --disable-build-doc --disable-starpufft-examples  --disable-fortran --with-perf-model-dir=$PWD/starpu_install 
                else
                        CFLAGS=-w ./configure --disable-starpufft --disable-cuda --disable-opencl --prefix=$PWD/starpu_install --enable-blas-lib=none --disable-starpu-top --disable-starpufft --disable-build-doc --disable-starpufft-examples --disable-mpi  --disable-fortran --disable-glpk --with-perf-model-dir=$PWD/starpu_install 
                fi
        fi
        make || make VERBOSE=1 || { echo 'STARPU installation failed' ; exit 1; }
        make install
        export STARPUROOT=$PWD
        export PKG_CONFIG_PATH=$STARPUROOT/starpu_install/lib/pkgconfig:$PKG_CONFIG_PATH
        export LD_LIBRARY_PATH=$STARPUROOT/starpu_install/lib:$LD_LIBRARY_PATH
        export CPATH=$STARPUROOT/starpu_install/include/starpu/1.2:$CPATH
        echo 'export PKG_CONFIG_PATH='$STARPUROOT'/starpu_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
        echo 'export LD_LIBRARY_PATH='$STARPUROOT'/starpu_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
        echo 'export CPATH='$STARPUROOT'/starpu_install/include/starpu/1.2:$CPATH' >> $SETUP_DIR/pkg_config.sh
    else
        echo "####################"
        echo "StarPU NOT FOUND"
        echo "Please download it from: http://starpu.gforge.inria.fr/files/"
        echo "After installing it, set the proper PKG_CONFIG_PATH variable"
        echo ""
        err=1
fi
fi

module load cmake/3.13.4/gnu-6.4.0
# CHAMELEON
if pkg-config --exists chameleon
then
    _LOCATION=`pkg-config --variable=prefix chameleon`
    echo "CHAMELEON FOUND in [$_LOCATION]"
else
    if [ "$BUILD_DEPENDENCIES" == "true" ]
    then
        echo "Building CHAMELEON..."
        cd $CHAMELEONDIR
        mkdir -p build/install_dir && cd build
        rm -rf ./CMake*
        cmake -DCHAMELEON_USE_MPI=$MPIVALUE -DCMAKE_BUILD_TYPE="Release" -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w" -DCHAMELEON_USE_CUDA=$CUDAVALUE -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$PWD/install_dir  -DMPI_C_COMPILER=`which mpicc` ..
        make || make VERBOSE=1 || { echo 'CHAMELEON installation failed' ; exit 1; }
        make install
    export PKG_CONFIG_PATH=$CHAMELEONDIR/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$CHAMELEONDIR/build/install_dir/lib:$LD_LIBRARY_PATH
    export CPATH=$CHAMELEONDIR/build/install_dir/include/coreblas:$CPATH
        export CPATH=$CHAMELEONDIR/build/install_dir/include:$CPATH
    echo 'export PKG_CONFIG_PATH='$CHAMELEONDIR'/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH'  >>  $SETUP_DIR/pkg_config.sh
    echo 'export LD_LIBRARY_PATH='$CHAMELEONDIR'/build/install_dir/lib:$LD_LIBRARY_PATH'  >>  $SETUP_DIR/pkg_config.sh
    echo 'export CPATH='$CHAMELEONDIR'/build/install_dir/include/coreblas:$CPATH'  >>  $SETUP_DIR/pkg_config.sh
        echo 'export CPATH='$CHAMELEONDIR'/build/install_dir/include:$CPATH'  >>  $SETUP_DIR/pkg_config.sh
    else
        echo "####################"
        echo "CHAMELEON NOT FOUND"
        echo "Please download it from: https://gitlab.inria.fr/solverstack/chameleon.git"
        echo "After installing it, set the proper PKG_CONFIG_PATH variable"
        echo ""
        err=1
    fi
fi

# starsh
if pkg-config --exists --atleast-version=0.1.1 starsh
then
    _LOCATION=`pkg-config --variable=prefix starsh`
    echo "starsh FOUND in [$_LOCATION]"
else
    if [ "$BUILD_DEPENDENCIES" == "true" ]
    then
        echo "Building STARS-H..."
        cd $STARSHDIR
    rm -rf build
        mkdir -p build/install_dir && cd build
        rm -rf ./CMake*
        cmake -DCMAKE_C_FLAGS=-fPIC -DCMAKE_BUILD_TYPE="Release" -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w" -DOPENMP=OFF -DSTARPU=OFF  -DEXAMPLES=OFF -DTESTING=OFF -DMPI=$MPIVALUE -DCMAKE_INSTALL_PREFIX=$PWD/install_dir   -DGSL=ON ..
        make || make VERBOSE=1 || { echo 'STARS-H installation failed' ; exit 1; }
        make install
    export PKG_CONFIG_PATH=$STARSHDIR/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$STARSHDIR/build/install_dir/lib:$LD_LIBRARY_PATH

    echo 'export PKG_CONFIG_PATH='$STARSHDIR'/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
    echo 'export PKG_CONFIG_PATH='$STARSHDIR'/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH' >>  $SETUP_DIR/pkg_config.sh
    else
        echo "####################"
        echo "starsh NOT FOUND"
        echo "Please download it from: https://github.com/ecrc/stars-h"
        echo "After installing it, set the proper PKG_CONFIG_PATH variable"
        echo ""
        err=1
    fi
fi

module unload cmake

module load cmake/3.9.4/gnu-6.4.0
export PATH=$PATH:/sw/csgv/openmpi/4.0.1/el7.6_gnu6.4.0_cuda10.1.105/bin
export CPATH=$CPATH:/sw/csgv/openmpi/4.0.1/el7.6_gnu6.4.0_cuda10.1.105/bin
# hicma
if pkg-config --exists  hicma
then
    _LOCATION=`pkg-config --variable=prefix hicma`
    echo "hicma FOUND in [$_LOCATION]"
else
    if [ "$BUILD_DEPENDENCIES" == "true" ]
    then
        echo "Building HiCMA..."
        cd $HICMADIR
    rm -rf build
        mkdir -p build/install_dir && cd build
        rm -rf ./CMake*
        cmake -DHICMA_USE_MPI=$MPIVALUE -DCMAKE_BUILD_TYPE="Release" -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w" -DBUILD_SHARED_LIBS=ON -DHICMA_ENABLE_TESTING=OFF -DHICMA_ENABLE_TIMING=OFF -DCMAKE_INSTALL_PREFIX=$PWD/install_dir   -DMPI_C_COMPILER=`which mpicc` .. 
        make || make VERBOSE=1 || { echo 'HICMA installation failed' ; exit 1; }
        make install
    export PKG_CONFIG_PATH=$HICMADIR/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$HICMADIR/build/install_dir/lib:$LD_LIBRARY_PATH


    echo 'export PKG_CONFIG_PATH='$HICMADIR'/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH' >>  $SETUP_DIR/pkg_config.sh
    echo 'export LD_LIBRARY_PATH='$HICMADIR'/build/install_dir/lib:$LD_LIBRARY_PATH' >>  $SETUP_DIR/pkg_config.sh
    else
        echo "####################"
        echo "hicma NOT FOUND"
        echo "Please download it from: https://github.com/ecrc/hicma"
        echo "After installing it, set the proper PKG_CONFIG_PATH variable"
        echo ""
        err=1
    fi
fi

# exageostat
if pkg-config --exists  exageostat
then
    _LOCATION=`pkg-config --variable=prefix exageostat`
    echo "exageostat FOUND in [$_LOCATION]"
else
    if [ "$BUILD_DEPENDENCIES" == "true" ]
    then
        echo "Building exageostat..."
        cd $EXAGEOSTATDEVDIR
    rm -rf build
        mkdir -p build/install_dir && cd build
        rm -rf ./CMake*
     cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w" -DBUILD_SHARED_LIBS=ON -DEXAGEOSTAT_EXAMPLES=ON   -DEXAGEOSTAT_USE_MPI=$MPIVALUE -DEXAGEOSTAT_USE_HICMA=ON -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DEXAGEOSTAT_USE_CUDA=$CUDAVALUE ..

        make || make VERBOSE=1 || { echo 'HICMA installation failed' ; exit 1; }
        make install
    else
        echo "####################"
        echo "hicma NOT FOUND"
        echo "Please download it from: https://github.com/ecrc/hicma"
        echo "After installing it, set the proper PKG_CONFIG_PATH variable"
        echo ""
        err=1
    fi
fi