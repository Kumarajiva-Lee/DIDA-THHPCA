#!/bin/sh
### global settings ######################
DIDAROOT=$PWD/..
#. $DIDAROOT/SCRIPTS/setenv.sh

CMP=$1

build_dir=build
threads=24

if [[ $1 == 'intel' ]]; then
  export CC=mpiicc CXX=sw9g++ FC=mpiifort
fi

if [[ $1 == 'sw' ]]; then
  export CC=mpicc CXX=mpicxx FC=mpif90
fi

### compile ######################################
if [ $# -eq 0 ]; then
  if [[ -d ${DIDAROOT}/$build_dir ]]; then
    cd ${DIDAROOT}/$build_dir
    make -j $threads
  else
    echo "Usage: "
    echo "compile.sh COMPILER(intel sw)"
    exit
  fi

  if [ -e dida.exe ]; then
    echo '     'dida.exe created!
  else
    echo '     'dida.exe NOT created!
  fi
  sleep 1

  #################################################
else
  ### compile ######################################
  echo ' '
  echo ======================================================================
  echo  compile process for whole dida system start.
  echo ======================================================================

  if [[ -d ${DIDAROOT}/$build_dir ]]; then
    rm -r ${DIDAROOT}/$build_dir
  fi
  mkdir ${DIDAROOT}/$build_dir
  cd ${DIDAROOT}/$build_dir
  cmake ..
  make -j $threads

  if [ -e dida.exe ]; then
    echo '     'dida.exe created!
  else
    echo '     'dida.exe NOT created!
  fi
  sleep 1

  echo ' '
  echo ======================================================================
  echo  compile process for whole dida system complete.
  echo ======================================================================
fi

