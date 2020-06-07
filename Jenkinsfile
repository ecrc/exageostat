pipeline {
/*
 * Defining where to run
 */
//// Any:
// agent any
//// By agent label:
//    agent { label 'sandybridge' }

    agent { label 'jenkinsfile' }
    triggers {
        pollSCM('H/15 * * * *')
    }
    environment {
        XX="gcc"
    }

    options {
        disableConcurrentBuilds()
        buildDiscarder(logRotator(numToKeepStr: '100'))
        timestamps()
    }

    stages {
        stage ('build') {
            steps {
                sh '''#!/bin/bash -le
module purge
module load ecrc-extras
module load mkl/2018-update-1
module load gcc/5.5.0
module load cmake/3.11.1
module load gsl/2.4-gcc-5.5.0
module load nlopt/2.4.2-gcc-5.5.0
module load hwloc/1.11.8-gcc-5.5.0
module load starpu/1.2.3-gcc-5.5.0-mkl-openmpi-3.0.0
module load hdf5/1.10.1-gcc-5.5.0
module load netcdf/4.5.0-gcc-5.5.0


module list
set -x

export EXAGEOSTATDEVDIR=$PWD
export HICMADIR=$EXAGEOSTATDEVDIR/hicma
export CHAMELEONDIR=$HICMADIR/chameleon
export STARSHDIR=$EXAGEOSTATDEVDIR/stars-h

## STARS-H
cd $STARSHDIR
rm -rf build
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DMPI=OFF -DOPENMP=OFF -DSTARPU=OFF
make clean
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

## CHAMELEON
cd $CHAMELEONDIR
# Update submodules
git submodule update --init --recursive
#install Chameleon
cd $CHAMELEONDIR
rm -rf build
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCHAMELEON_USE_MPI=OFF  -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DCHAMELEON_ENABLE_EXAMPLE=OFF -DCHAMELEON_ENABLE_TESTING=OFF -DCHAMELEON_ENABLE_TIMING=OFF
make -j # CHAMELEON parallel build seems to be fixed
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

## HICMA
cd $HICMADIR
rm -rf build
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DHICMA_USE_MPI=OFF -DHICMA_ENABLE_TESTING=OFF -DHICMA_ENABLE_TIMING=OFF
make clean
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

# EXAGEOSTAT
cd $EXAGEOSTATDEVDIR
rm -rf build
mkdir -p build
cd build
cmake .. \
    -DCMAKE_INSTALL_PREFIX=$PWD/installdir \
    -DEXAGEOSTAT_SCHED_STARPU=ON \
    -DEXAGEOSTAT_USE_MPI=OFF \
    -DEXAGEOSTAT_PACKAGE=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DEXAGEOSTAT_USE_STARSH=ON \
    -DEXAGEOSTAT_USE_HICMA=ON \
    -DEXAGEOSTAT_USE_NETCDF=ON

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
module load mkl/2018-update-1
module load gcc/5.5.0
module load cmake/3.11.1
module load gsl/2.4-gcc-5.5.0
module load nlopt/2.4.2-gcc-5.5.0
module load hwloc/1.11.8-gcc-5.5.0
module load starpu/1.2.3-gcc-5.5.0-mkl-openmpi-3.0.0
module load hdf5/1.10.1-gcc-5.5.0
module load netcdf/4.5.0-gcc-5.5.0


module list
set -x

cd build
rm -rf Testing
ctest --no-compress-output -T Test

'''
                step([$class: 'XUnitBuilder',
                     thresholds: [[$class: 'FailedThreshold', unstableThreshold: '0']],
                     tools: [[$class: 'CTestType', pattern: 'build/Testing/**/Test.xml']]])
            }
        }
        stage ('docs') {
            steps {
                sh "cd $WORKSPACE/build && make docs"
                sh '''#!/bin/bash -ex
                      cd $WORKSPACE
                      rm -rf cppcheckhtml
                      cppcheck --enable=all --xml --xml-version=2 exageostat_approx/ exageostat_exact/ examples/ misc/ src/ r-wrappers/ -I include/   2> cppcheck.xml
                      cppcheck-htmlreport --source-encoding="iso8859-1" --title="ExaGeoStat" --source-dir=. --report-dir=cppcheckhtml --file=cppcheck.xml
'''
                publishHTML( target: [allowMissing: false, alwaysLinkToLastBuild: false, keepAll: false, reportDir: 'cppcheckhtml', reportFiles: 'index.html', reportName: 'CppCheckReport', reportTitles: ''] )
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

    // Post build actions
    post {
        //always {
        //}
        //success {
        //}
        //unstable {
        //}
        //failure {
        //}
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


