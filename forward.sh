#!/bin/bash

echo "******************* running example for the z coponent: `date` ********************"
currentdir=`pwd`
rm -rf adt DATA SEM x_result nohup.out z_result
mkdir -p OUTPUT_FILES
mkdir -p DATA
mkdir -p SEM

# sets up local DATA/ directory
cd DATA/
ln -s ../Par_file_zh Par_file
ln -s ../SOURCE_zh SOURCE
ln -s ../interfaces_zh.dat interfaces.dat
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

cd $currentdir

# links executables
rm -f xmeshfem2D xspecfem2D xsmooth_sem
ln -s ../../bin/xmeshfem2D
ln -s ../../bin/xspecfem2D
ln -s ../../bin/xsmooth_sem

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running mesher..."
  echo
  ./xmeshfem2D
else
  # This is a MPI simulation
  echo
  echo "running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC ./xmeshfem2D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running solver..."
  echo
  ./xspecfem2D
else
  # This is a MPI simulation
  echo
  echo "running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./xspecfem2D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES
echo ">>>>>>>>>>>> the forward problem with x is done !!!"

rm OUTPUT_FILES/pml*
gfortran adj_seismogram_AH.f90 constants.f90 fft.f90 gauss.f90 m_hilbert_transform.f90 window.f90 -o adt
./adt
