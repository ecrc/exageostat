#!/usr/bin/env bash

# Installation script for StarPU

#VARIABLES
CUDA_VALUE="OFF"
MPI_VALUE="OFF"

print_usage() {
	echo "usage: $0 --prefix /path/to/install --setup /path/to/download/tar/file [--cuda ON|OFF] [--mpi ON|OFF]"
}

while [ -n "$1"  ]
do
	case "$1" in
		--cuda)
		  shift
			CUDA_VALUE=$1
			shift
			;;
		--mpi)
		  shift
			MPI_VALUE=$1
			shift
			;;
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

echo "Installing StarPu..."
cd "$SETUP_DIR" || exit
wget https://files.inria.fr/starpu/starpu-1.3.9/starpu-1.3.9.tar.gz -O - | tar -zx
echo "Building StarPu... with CUDA $CUDA_VALUE AND MPI $MPI_VALUE"
cd starpu* || exit
if [ "$CUDA_VALUE" == "ON" ]; then
  if [ "$MPI_VALUE" == "ON" ]; then
    ./configure --disable-starpufft --enable-cuda --disable-opencl --prefix="$PREFIX"   --disable-starpu-top --disable-starpufft --disable-build-doc --disable-starpufft-examples   --disable-fortran --with-perf-model-dir="$SETUP_DIR"  --disable-fstack-protector-all --disable-gcc-extensions
  else
    ./configure --disable-starpufft --enable-cuda  --disable-mpi --disable-opencl --prefix="$PREFIX"   --disable-starpu-top --disable-starpufft --disable-build-doc --disable-starpufft-examples   --disable-fortran --with-perf-model-dir="$SETUP_DIR" --disable-fstack-protector-all --disable-gcc-extensions
  fi
else
  if [ "$MPI_VALUE" == "ON" ]; then
    ./configure --disable-starpufft --disable-cuda --disable-opencl --prefix="$PREFIX"  --disable-starpu-top --disable-starpufft --disable-build-doc --disable-starpufft-examples  --disable-fortran --with-perf-model-dir="$SETUP_DIR" --disable-fstack-protector-all --disable-gcc-extensions
  else
    ./configure --disable-starpufft --disable-cuda --disable-mpi --disable-opencl --prefix="$PREFIX" --disable-starpu-top --disable-starpufft --disable-build-doc --disable-starpufft-examples   --disable-fortran --disable-glpk --with-perf-model-dir="$SETUP_DIR" --disable-fstack-protector-all --disable-gcc-extensions
  fi
fi
make clean && make -j || { echo 'STARPU installation failed' ; exit 1; }
make -j install

export PKG_CONFIG_PATH=$PREFIX/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH