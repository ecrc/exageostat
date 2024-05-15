pipeline {
    agent none
    triggers {
        pollSCM('H/15 * * * *')
    }
    options {
        disableConcurrentBuilds()
        buildDiscarder(logRotator(numToKeepStr: '100'))
        timestamps()
    }

    stages {
        stage ('with-HiCMA'){
            agent { label 'jenkinsfile' }
            stages{
                stage ('build') {
                    steps {
                        sh '''#!/bin/bash -le
                            module purge
                            module load ecrc-extras
                            module load mkl/2020.0.166
                            module load gcc/10.2.0
                            module load cmake/3.21.2
                            module load hwloc/2.4.0-gcc-10.2.0
                            module load openmpi/4.1.0-gcc-10.2.0
                            module load starpu/1.3.9-gcc-10.2.0-mkl-openmpi-4.1.0
                            module load gsl/2.6-gcc-10.2.0
                            module load nlopt/2.7.0-gcc-10.2.0
                            module load hdf5/1.12.0-gcc-10.2.0
                            module load netcdf/4.7.4-gcc-10.2.0
                            module load doxygen/1.8.20

                            module list
                            set -x

                            export EXAGEOSTATDEVDIR=$PWD
                            export HICMADIR=$EXAGEOSTATDEVDIR/hicma
                            export CHAMELEONDIR=$EXAGEOSTATDEVDIR/chameleon
                            export STARSHDIR=$HICMADIR/stars-h
                            export HCOREDIR=$HICMADIR/hcore
                            export HICMAINSTALLDIR=$HICMADIR/dependencies-prefix

                            ## CHAMELEON
                            # Update submodules
                            git submodule update --init --recursive
                            cd $CHAMELEONDIR
                            git checkout release-1.1.0

                            #install Chameleon
                            rm -rf build
                            mkdir -p build/installdir
                            cd build
                            cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DCMAKE_C_FLAGS=-fPIC -DCHAMELEON_USE_MPI=OFF -DCMAKE_BUILD_TYPE="Release" \
                            -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w" -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_ENABLE_EXAMPLE=OFF \
                            -DCHAMELEON_ENABLE_TESTING=OFF -DCHAMELEON_ENABLE_TIMING=OFF -DBUILD_SHARED_LIBS=OFF

                            make clean
                            make -j # CHAMELEON parallel build seems to be fixed
                            make install
                            export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

                            ## HICMA
                            cd $HICMADIR
                            git submodule update --init --recursive
                            mkdir -p $HICMAINSTALLDIR
                            rm -rf $HICMAINSTALLDIR/*

                                # STARS-H
                                cd $STARSHDIR
                                rm -rf build
                                mkdir -p build
                                cd build
                                cmake .. -DCMAKE_INSTALL_PREFIX=$HICMAINSTALLDIR -DMPI=OFF -DOPENMP=OFF -DSTARPU=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_C_FLAGS=-fPIC -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w"
                                make clean
                                make -j
                                make install

                                ## HCORE
                                cd $HCOREDIR
                                rm -rf build
                                mkdir -p build
                                cd build
                                cmake .. -DCMAKE_INSTALL_PREFIX=$HICMAINSTALLDIR -DCMAKE_C_FLAGS=-fPIC -DBUILD_SHARED_LIBS=ON -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w"
                                make clean
                                make -j
                                make install

                            export PKG_CONFIG_PATH=$HICMAINSTALLDIR/lib/pkgconfig:$PKG_CONFIG_PATH
                            cd $HICMADIR
                            rm -rf build
                            mkdir -p build/installdir
                            cd build
                            cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DCMAKE_C_FLAGS=-fPIC -DHICMA_USE_MPI="$MPI_VALUE" -DCMAKE_BUILD_TYPE="Release" \
                            -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w" -DBUILD_SHARED_LIBS=ON -DCMAKE_C_FLAGS="-fcommon"
                            make clean
                            make -j
                            make install
                            export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

                            # EXAGEOSTAT
                            cd $EXAGEOSTATDEVDIR
                            rm -rf build
                            mkdir -p build/installdir
                            cd build
                            CFLAGS="-fcommon" cmake .. \
                                -DCMAKE_INSTALL_PREFIX=$PWD/installdir \
                                -DEXAGEOSTAT_SCHED_STARPU=ON \
                                -DEXAGEOSTAT_USE_MPI=OFF \
                                -DEXAGEOSTAT_PACKAGE=ON \
                                -DCMAKE_BUILD_TYPE=Release \
                                -DEXAGEOSTAT_USE_HICMA=ON \
                                -DEXAGEOSTAT_USE_NETCDF=OFF \
                                -DEXAGEOSTAT_USE_CHAMELEON=ON \
                                -DEXAGEOSTAT_INSTALL_DEPS=OFF \
                                -DMPI_VALUE=OFF \
                                -DEXAGEOSTAT_USE_CUDA=OFF \
                                -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w" \
                                -DBUILD_SHARED_LIBS=ON
                            make clean
                            make -j || make VERBOSE=1
                            make install
                            '''
                    }
                }
                stage ('test') {
                    steps {
                        sh '''#!/bin/bash -le
                            module purge
                            module load ecrc-extras
                            module load mkl/2020.0.166
                            module load gcc/10.2.0
                            module load cmake/3.21.2
                            module load hwloc/2.4.0-gcc-10.2.0
                            module load openmpi/4.1.0-gcc-10.2.0
                            module load starpu/1.3.9-gcc-10.2.0-mkl-openmpi-4.1.0
                            module load gsl/2.6-gcc-10.2.0
                            module load nlopt/2.7.0-gcc-10.2.0
                            module load hdf5/1.12.0-gcc-10.2.0
                            module load netcdf/4.7.4-gcc-10.2.0

                            module list
                            set -x

                            cd build
                            rm -rf Testing
                            ctest --no-compress-output -T Test
                            '''
                    }
                }
                stage ('docs') {
                    steps {
                        sh "cd $WORKSPACE/build && make docs"
                        publishHTML( target: [allowMissing: false, alwaysLinkToLastBuild: false, keepAll: false, reportDir: 'build/docs/build/html', reportFiles: 'index.html', reportName: 'Doxygen Documentation', reportTitles: ''] )
                    }
                }
                stage ('packages') {
                    steps {
                        sh "cd $WORKSPACE/build && make package"
                        archiveArtifacts allowEmptyArchive: true, artifacts: 'build/exageostat-*gz'
                    }
                }
            }
        }

        stage ('Without-HiCMA'){
            agent { label 'jenkinsfile' }
            stages {
                stage ('build') {
                    steps {
                        sh '''#!/bin/bash -le
                            module purge
                            module load ecrc-extras
                            module load mkl/2020.0.166
                            module load gcc/10.2.0
                            module load cmake/3.21.2
                            module load hwloc/2.4.0-gcc-10.2.0
                            module load openmpi/4.1.0-gcc-10.2.0
                            module load starpu/1.3.9-gcc-10.2.0-mkl-openmpi-4.1.0
                            module load gsl/2.6-gcc-10.2.0
                            module load nlopt/2.7.0-gcc-10.2.0
                            module load hdf5/1.12.0-gcc-10.2.0
                            module load netcdf/4.7.4-gcc-10.2.0
                            module load doxygen/1.8.20

                            module list
                            set -x

                            export EXAGEOSTATDEVDIR=$PWD
                            export CHAMELEONDIR=$EXAGEOSTATDEVDIR/chameleon

                            # Update submodules
                            git submodule update --init --recursive

                            ## CHAMELEON
                            cd $CHAMELEONDIR
                            #install Chameleon
                            rm -rf build
                            mkdir -p build/installdir
                            cd build
                            cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DCMAKE_C_FLAGS=-fPIC -DCHAMELEON_USE_MPI=OFF -DCMAKE_BUILD_TYPE="Release" \
                            -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w" -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_ENABLE_EXAMPLE=OFF \
                            -DCHAMELEON_ENABLE_TESTING=OFF -DCHAMELEON_ENABLE_TIMING=OFF -DBUILD_SHARED_LIBS=OFF

                            make clean
                            make -j
                            make install
                            export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

                            # EXAGEOSTAT
                            cd $EXAGEOSTATDEVDIR
                            rm -rf build
                            mkdir -p build/installdir
                            cd build
                            CFLAGS="-fcommon" cmake .. \
                                -DCMAKE_INSTALL_PREFIX=$PWD/installdir \
                                -DEXAGEOSTAT_SCHED_STARPU=ON \
                                -DEXAGEOSTAT_USE_MPI=OFF \
                                -DEXAGEOSTAT_PACKAGE=ON \
                                -DCMAKE_BUILD_TYPE=Release \
                                -DEXAGEOSTAT_USE_HICMA=OFF \
                                -DEXAGEOSTAT_USE_NETCDF=ON \
                                -DEXAGEOSTAT_USE_CHAMELEON=ON \
                                -DEXAGEOSTAT_INSTALL_DEPS=OFF \
                                -DMPI_VALUE=OFF \
                                -DCUDA_VALUE=OFF \
                                -DCMAKE_C_FLAGS_RELEASE="-O3 -Ofast -w" \
                                -DBUILD_SHARED_LIBS=ON
                            make clean
                            make -j || make VERBOSE=1
                            make install
                            '''
                    }
                }
                stage ('test') {
                    steps {
                        sh '''#!/bin/bash -le
                            module purge
                            module load ecrc-extras
                            module load mkl/2020.0.166
                            module load gcc/10.2.0
                            module load cmake/3.21.2
                            module load hwloc/2.4.0-gcc-10.2.0
                            module load openmpi/4.1.0-gcc-10.2.0
                            module load starpu/1.3.9-gcc-10.2.0-mkl-openmpi-4.1.0
                            module load gsl/2.6-gcc-10.2.0
                            module load nlopt/2.7.0-gcc-10.2.0
                            module load hdf5/1.12.0-gcc-10.2.0
                            module load netcdf/4.7.4-gcc-10.2.0

                            module list
                            set -x

                            cd build
                            rm -rf Testing
                            ctest --no-compress-output -T Test

                            '''
                    }
                }
            }
        }

    }

    // Post build actions
    post {
        unstable {
            emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build is UNSTABLE", recipientProviders: [culprits(),requestor()]
        }
        fixed {
            emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build FIXED!", recipientProviders: [culprits(),requestor()]
        }
        failure {
            emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build FAILED", recipientProviders: [culprits(),requestor()]
        }
    }
}