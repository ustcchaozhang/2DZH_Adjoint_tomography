#!/bin/bash
echo "******************* running example for the x coponent: `date` ********************"
currentdir=`pwd`
rm -rf nohup.out a.out SEM SEM_X SEM_Z DATA_X DATA_Z x_result z_result OUTPUT_FILES
mkdir -p OUTPUT_FILES
mkdir -p DATA
mkdir -p SEM

sed -i "47s:.*:adj_comp = 1:g" adj_seismogram_AH.f90
sed -i "13s:.*:SIMULATION_TYPE          = 1:g" Par_file_zh
sed -i "17s:.*:SAVE_FORWARD             =.true.:g" Par_file_zh

# sets up local DATA/ directory
cd DATA/
ln -s ../Par_file_zh Par_file
ln -s ../SOURCE_zh SOURCE
ln -s ../interfaces_zh.dat interfaces.dat
cd ../

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
  echo "running mesher..."
  ./xmeshfem2D > OUTPUT_FILES/output_mesher.txt
else
  # This is a MPI simulation
  echo "running mesher on $NPROC processors..."
  mpirun -np $NPROC ./xmeshfem2D > OUTPUT_FILES/output_mesher.txt
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo "running solver..."
  ./xspecfem2D > OUTPUT_FILES/output_solver.txt
else
  # This is a MPI simulation
  echo "running solver on $NPROC processors..."
  mpirun -np $NPROC ./xspecfem2D > OUTPUT_FILES/output_solver.txt
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES
echo ">>>>>>>>>>>> the forward problem with x is done !!!"

gfortran adj_seismogram_AH.f90 constants.f90 fft.f90 gauss.f90 m_hilbert_transform.f90 window.f90 -o adt
./adt >OUTPUT_FILES/output_adjoint.txt

sed -i "13s:.*:SIMULATION_TYPE          = 3:g" Par_file_zh
sed -i "17s:.*:SAVE_FORWARD             =.false.:g" Par_file_zh
 
# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo "running mesher..."
  ./xmeshfem2D > OUTPUT_FILES/output_mesher_kernel.txt
else
  # This is a MPI simulation
  echo "running mesher on $NPROC processors..."
  mpirun -np $NPROC ./xmeshfem2D > OUTPUT_FILES/output_mesher_kernel.txt
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo "running solver..."
  ./xspecfem2D > OUTPUT_FILES/output_solver_kernel.txt
else
  # This is a MPI simulation
  echo "running solver on $NPROC processors..."
  mpirun -np $NPROC ./xspecfem2D > OUTPUT_FILES/output_solver_kernel.txt
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo ">>>>>>>>>>>> the kernel problem with x is done !!!"

sed -i "13s:.*:SIMULATION_TYPE          = 1:g" Par_file_zh
sed -i "17s:.*:SAVE_FORWARD             =.true.:g" Par_file_zh

mv ./DATA/*.bin ./OUTPUT_FILES
echo "running smooth on $NPROC processors..."
mpirun -np $NPROC ./xsmooth_sem 12000 6000 0 beta_kernel OUTPUT_FILES OUTPUT_FILES false > OUTPUT_FILES/output_smooth.txt

mv OUTPUT_FILES x_result
mv DATA DATA_X
mv SEM SEM_X
echo ">>>>>>>>>>>> x result saved in x_reslut !!!"

echo "******************** running example for the z coponent: `date` ********************"

sed -i "47s:.*:adj_comp = 3:g" adj_seismogram_AH.f90

mkdir -p OUTPUT_FILES
mkdir -p DATA
mkdir -p SEM

# sets up local DATA/ directory
cd DATA/
ln -s ../Par_file_zh Par_file
ln -s ../SOURCE_zh SOURCE
ln -s ../interfaces_zh.dat interfaces.dat
cd ../

cd $currentdir

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# links executables
rm -f xmeshfem2D xspecfem2D xsmooth_sem
ln -s ../../bin/xmeshfem2D
ln -s ../../bin/xspecfem2D
ln -s ../../bin/xsmooth_sem

# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo "running mesher..."
  ./xmeshfem2D > OUTPUT_FILES/output_mesher.txt
else
  # This is a MPI simulation
  echo "running mesher on $NPROC processors..."
  mpirun -np $NPROC ./xmeshfem2D > OUTPUT_FILES/output_mesher.txt
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo "running solver..."
  ./xspecfem2D > OUTPUT_FILES/output_solver.txt
else
  # This is a MPI simulation
  echo "running solver on $NPROC processors..."
  mpirun -np $NPROC ./xspecfem2D > OUTPUT_FILES/output_solver.txt
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES
echo ">>>>>>>>>>>> the forward problem with z is done !!!"

gfortran adj_seismogram_AH.f90 constants.f90 fft.f90 gauss.f90 m_hilbert_transform.f90 window.f90 -o adt
./adt > OUTPUT_FILES/output_adjoint.txt

sed -i "13s:.*:SIMULATION_TYPE          = 3:g" Par_file_zh
sed -i "17s:.*:SAVE_FORWARD             =.false.:g" Par_file_zh

# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo "running mesher..."
  ./xmeshfem2D > OUTPUT_FILES/output_mesher_kernel.txt
else
  # This is a MPI simulation
  echo "running mesher on $NPROC processors..."
  mpirun -np $NPROC ./xmeshfem2D > OUTPUT_FILES/output_mesher_kernel.txt
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo "running solver..."
  ./xspecfem2D > OUTPUT_FILES/output_solver_kernel.txt
else
  # This is a MPI simulation
  echo "running solver on $NPROC processors..."
  mpirun -np $NPROC ./xspecfem2D > OUTPUT_FILES/output_solver_kernel.txt
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo ">>>>>>>>>>>> the kernel problem with z is done !!!"

sed -i "13s:.*:SIMULATION_TYPE          = 1:g" Par_file_zh
sed -i "17s:.*:SAVE_FORWARD             =.true.:g" Par_file_zh

mv ./DATA/*.bin ./OUTPUT_FILES
echo "running smooth on $NPROC processors..."
mpirun -np $NPROC ./xsmooth_sem 12000 6000 0 beta_kernel OUTPUT_FILES OUTPUT_FILES false > OUTPUT_FILES/output_smooth.txt

mv OUTPUT_FILES z_result
mv DATA DATA_Z
mv SEM SEM_Z
echo ">>>>>>>>>>>> z result saved in z_reslut !!!"
sed -i "47s:.*:adj_comp = 1:g" adj_seismogram_AH.f90

echo "******************* running example for the xz all done: `date` ********************"

rm x_result/pml*
rm z_result/pml*

python zh_plot.py x_result/ beta_kernel_smooth $NPROC z_result/
