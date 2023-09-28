#!/bin/bash
# sigmax_array=( 000 002 004 008 016 )
# sigmat_array=( 000 006 030 180 )
sigmax_array=( 000 )
sigmat_array=( 000 )
pathCases="../cases/"
pathOpenBC="${pathCases}openBC/"
cd $pathOpenBC
for sigmat in ${sigmat_array[@]}; do
	for sigmax in ${sigmax_array[@]}; do
		experiment="x${sigmax}y${sigmax}z000t${sigmat}"
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
done
cd ../../scripts/