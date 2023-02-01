#VARIABLES
MPI_VALUE="OFF"
CUDA_VALUE="OFF"

print_usage() {
	echo "usage: $0 --prefix /path/to/install --setup /path/to/download/tar/file [--mpi ON|OFF]"
}

echo "Reading Arguments"
while [ -n "$1"  ]
do
	case "$1" in
		--prefix)
		  shift
			PREFIX=$1
			shift
			;;
    --setup)
      shift
      SETUP_DIR=$1
      shift
      ;;
    --build)
      shift
      BUILD_DIR=$1
      shift
      ;;
    --mpi)
      shift
      MPI_VALUE=$1
      shift
      ;;
    --cuda)
      shift
      CUDA_VALUE=$1
      shift
      ;;
		--help|-h)
			print_usage
			exit 0
			;;
		*)
			print_usage
			exit 1
			;;
	esac
done

if [ "$MPI_VALUE" = "ON" ]
   	then
   	  MPI_WORKS="-DMPI_C_WORKS=ON"   #Hicma Sets Linker flags that interfere with CMake MPI build tests. This flag needs
   	                                 # to be set in advance to ensure MPI checks pass. MPI is already found and configured
   	                                 # at an earlier stage
fi

echo "Installing HiCMA dependencies..."
cd "$SETUP_DIR" || exit 1
git clone https://github.com/ecrc/hicma hicma
cd hicma || exit 1;
git checkout fa8596b5d3aa8e5b7d5c06cd8db3cecc32f70d17 #V1.0.0
git submodule update --init --recursive

echo "Building Stars-H"
cd stars-h || exit 1;
rm -rf build && mkdir -p build
cd build && cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DMPI=OFF -DOPENMP=OFF -DSTARPU=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_C_FLAGS=-fPIC -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w"
make clean && make -j || { echo 'Stars-H from HICMA installation failed' ; exit 1; }
make -j install

echo "Building HCORE"
cd ../../hcore || exit 0;
rm -rf build && mkdir -p build
cd build && cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_C_FLAGS=-fPIC -DBUILD_SHARED_LIBS=ON -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w"
make clean && make -j || { echo 'HCORE from HICMA installation failed' ; exit 1; }
make -j install

echo "Building HiCMA"
cd ../../ && rm -rf build && mkdir -p build
cd build || exit 1;
cp -r "$BUILD_DIR"/CMakeFiles .
cmake -DCMAKE_C_FLAGS=-fPIC -DHICMA_USE_MPI="$MPI_VALUE" -DCMAKE_BUILD_TYPE="Release" -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w" \
-DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_C_FLAGS="-fcommon" ..
make clean && make -j || { echo 'HICMA installation failed' ; exit 1; }
make -j install

export PKG_CONFIG_PATH=$PREFIX/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
