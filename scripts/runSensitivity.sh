#!/bin/bash
tauh_array=( 000 060 )
pbc_array=( 2 4 )
dint_array=( 00 05 )
pathCases="../cases/"
pathSensitivity="${pathCases}sensitivity/"
cd $pathSensitivity
# Create simulations for variation in tauh
for tauh in ${tauh_array[@]}; do
  experiment="tauh${tauh}"
  cd $experiment
  mpirun -np 8 ./dales4.4 2>&1 | tee output_sim.txt
  ./merge.sh 2>&1 | tee output_merge.txt 
  rm initd* 
  rm fielddump.00*
  rm field.001 
  rm flux*
  rm meancrossxz.00*
  rm moments.001 
  rm tmser1.001 
  rm tmsurf.001
  cd ..
done
# Create simulations for variation in pbc
for pbc in ${pbc_array[@]}; do
  experiment="pbc${pbc}"
  cd $experiment
  mpirun -np 8 ./dales4.4 2>&1 | tee output_sim.txt
  ./merge.sh 2>&1 | tee output_merge.txt 
  rm initd* 
  rm fielddump.00*
  rm field.001 
  rm flux*
  rm meancrossxz.00*
  rm moments.001 
  rm tmser1.001 
  rm tmsurf.001
  cd ..
done
# Create simulations for variation in dint
for dint in ${dint_array[@]}; do
  experiment="dint${dint}"
  cd $experiment
  mpirun -np 8 ./dales4.4 2>&1 | tee output_sim.txt
  ./merge.sh 2>&1 | tee output_merge.txt 
  rm initd* 
  rm fielddump.00*
  rm field.001 
  rm flux*
  rm meancrossxz.00*
  rm moments.001 
  rm tmser1.001 
  rm tmsurf.001
  cd ..
done
# Create simulation without buoyancy term at the top boundary
experiment="lbuoyoff"
cd $experiment
mpirun -np 8 ./dales4.4 2>&1 | tee output_sim.txt
./merge.sh 2>&1 | tee output_merge.txt 
rm initd* 
rm fielddump.00*
rm field.001 
rm flux*
rm meancrossxz.00*
rm moments.001 
rm tmser1.001 
rm tmsurf.001
cd ..
cd ../../scripts/