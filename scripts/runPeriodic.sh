#!/bin/bash
pathCases="../cases/"
pathPeriodic="${pathCases}periodic/"
cd $pathPeriodic
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
cd ../../scripts/

