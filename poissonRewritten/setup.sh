#!/bin/bash
mkdir release
cd release

module load cmake
module load intelcomp
module load openmpi/1.4.3-intel

CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make

#PBS -A freecycle 

cp ../job.sh ./job.sh
chmod 755 job.sh
