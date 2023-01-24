#!/usr/bin/env bash
# Installation script for GSL

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

echo "Installing GSL..."
cd "$SETUP_DIR" || exit
wget https://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz -O - | tar -zx
echo "Building GSL..."
cd gsl* || exit
./configure --prefix="$PREFIX"
make clean && make -j || { echo 'GSL installation failed' ; exit 1; }
make -j install

export PKG_CONFIG_PATH=$PREFIX/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH