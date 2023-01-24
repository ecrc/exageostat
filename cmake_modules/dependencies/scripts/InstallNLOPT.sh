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

if [ -z "$PREFIX" ] || [ -z "$SETUP_DIR" ]; then
  echo "Please specify Installation Path (Path where the library will be actually installed)
  AND Setup Path (Path where the untarred binaries will be downloaded)"
  exit 1;
fi

echo "Installing NLOPT..."
cd "$SETUP_DIR" || exit
wget https://github.com/stevengj/nlopt/archive/v2.7.0.tar.gz -O - | tar -zx
echo "Building NLOPT..."
cd nlopt* || exit
mkdir -p build
cd build || exit
rm -rf ./CMake*
cmake  -DCMAKE_C_FLAGS=-fcommon -DCMAKE_C_FLAGS=-fPIC -DCMAKE_BUILD_TYPE="Release" \
-DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w" -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX="$PREFIX"  ..
make clean && make -j || { echo 'NLOPT installation failed' ; exit 1; }
make -j install

export PKG_CONFIG_PATH=$PREFIX/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
export NLOPT_LIBRARY_DIRS=$PREFIX/lib:$NLOPT_LIBRARY_DIRS
export NLOPT_INCLUDE_DIRS=$PREFIX/include:NLOPT_INCLUDE_DIRS
