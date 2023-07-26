#!/bin/bash
# Create folders for simulations and copy initial profiles
sigmax_array=( 000 002 004 008 016 )
sigmat_array=( 000 006 030 180 )
pathSim="../cases/openBC_synturb/"
pathPer="../cases/periodic/"
mkdir $pathSim
for sigmat in ${sigmat_array[@]}; do
	for sigmax in ${sigmax_array[@]}; do
		experiment="x${sigmax}y${sigmax}z000t${sigmat}"
		# Skip simulation without smoothing as it's the same as the openBC simulation
    if [ $experiment = "x000y000z000t000" ]; then
      continue
    fi
		# create simulation directory
		mkdir "${pathSim}${experiment}"
		cp "${pathPer}prof.inp.000" "${pathSim}${experiment}/prof.inp.001"
		cp "${pathPer}lscale.inp.000" "${pathSim}${experiment}/lscale.inp.001"
		cp "namoptions.openBC_synturb" "${pathSim}${experiment}/namoptions" 
	done
done
