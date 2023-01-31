#!/usr/bin/env bash

# Installation script for HWLOC

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

echo "Installing HWLOC..."
cd "$SETUP_DIR" || exit
wget https://download.open-mpi.org/release/hwloc/v2.4/hwloc-2.4.0.tar.gz -O - | tar -zx
echo "Building HWLOC..."
cd hwloc* || exit
./configure --prefix="$PREFIX"
make clean && make -j || { echo 'HWLOC installation failed' ; exit 1; }
make -j install

export PKG_CONFIG_PATH=$PREFIX/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH