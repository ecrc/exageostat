Installation requires at least **CMake** of version 3.2.3. To build ExaGeoStat,
follow these instructions:

1.  Get  from git repository

        git clone git@github.com:ecrc/exageostat

    or

        git clone https://github.com/ecrc/exageostat

2.  Go into exageostat folder

        cd exageostat

3.  Get submodules

        git submodule update --init --recursive

4.  Create build directory and go there

        mkdir build && cd build

5.  Use CMake to get all the dependencies

        cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/  -DEXAGEOSTAT_SCHED_STARPU=ON   -DEXAGEOSTAT_USE_NETCDF=ON -DEXAGEOSTAT_USE_HICMA=ON


6.  Build EXAGEOSTAT

        make -j

8.  Build local documentation (optional)

        make docs

9.  Install EXAGEOSTAT

        make install

10. Add line

        export PKG_CONFIG_PATH=/path/to/install/lib/pkgconfig:$PKG_CONFIG_PATH

    to your .bashrc file.

Now you can use pkg-config executable to collect compiler and linker flags for
EXAGEOSTAT.


