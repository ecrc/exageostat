#!/usr/bin/env bash

# Installation script for NLOPT

print_usage() {
	echo "usage: $0 --prefix /path/to/install --setup /path/to/download/tar/file"
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

if [ -z "$PREFIX" ]; then
  echo "Please specify Installation Path (Path where the library will be actually installed)
  AND Setup Path (Path where the untarred binaries will be downloaded)"
  exit 1;
fi

echo "Installing hdf5 library for netcdf4..."
cd "$SETUP_DIR" || exit
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.gz -O - | tar -zx
echo "Building hdf5..."
cd hdf5*
./configure --prefix="$PREFIX" --with-zlib="$PREFIX" --enable-hl
make clean && make -j || { echo 'hdf5 installation failed' ; exit 1; }
make -j install
echo "hdf5 Done..."

export HDF5=$PREFIX
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH

echo "Installing NetCDF..."
cd "$SETUP_DIR" || exit
wget https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.7.4.tar.gz -O - | tar -zx
echo "Building NetCDF..."
cd netcdf* || exit

export CPPFLAGS=-I$PREFIX/include
export LDFLAGS=-L$PREFIX/lib

./configure --disable-dap --prefix="$PREFIX"
make clean && make -j || { echo 'Netcdf installation failed' ; exit 1; }
make -j install

export PKG_CONFIG_PATH=$PREFIX/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH