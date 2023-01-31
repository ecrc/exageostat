#VARIABLES
MPI_VALUE="OFF"
CUDA_VALUE="OFF"

print_usage() {
	echo "usage: $0 --prefix /path/to/install --setup /path/to/download/tar/file [--cuda ON|OFF] [--mpi ON|OFF]"
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

echo "Cloning Chameleon..."
cd "$SETUP_DIR" || exit 1
git clone https://gitlab.inria.fr/solverstack/chameleon.git chameleon
cd chameleon || exit 1;
git checkout release-1.1.0
git submodule update --init --recursive
echo "Building Chameleon..."
rm -rf build
mkdir -p build
cd build || exit 1
rm -rf ./CMake*

LDFLAGS="-L$PREFIX/lib" cmake -DCMAKE_C_FLAGS=-fPIC -DCHAMELEON_USE_MPI="$MPI_VALUE" -DCMAKE_BUILD_TYPE="Release" \
-DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w" -DCHAMELEON_USE_CUDA="$CUDA_VALUE" -DCHAMELEON_ENABLE_EXAMPLE=OFF \
-DCHAMELEON_ENABLE_TESTING=OFF -DCHAMELEON_ENABLE_TIMING=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="$PREFIX" ..

make clean && make -j || { echo 'CHAMELEON installation failed' ; exit 1; }
make -j install

export PKG_CONFIG_PATH=$PREFIX/installdir/lib/pkgconfig:$PKG_CONFIG_PATH