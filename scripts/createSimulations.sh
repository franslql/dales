#!/bin/bash
# Create folders for simulations and copy initial profiles
type="openBC"
sigmax_array=( 000 002 004 008 016 )
sigmat_array=( 000 006 030 180 )
pathSim="../cases/${type}/"
pathPer="../cases/periodic/"
pathDALES="/perm/nmfl/dales/build/src/"
mkdir $pathSim
for sigmat in ${sigmat_array[@]}; do
	for sigmax in ${sigmax_array[@]}; do
		experiment="x${sigmax}y${sigmax}z000t${sigmat}"
		# create simulation directory
		mkdir "${pathSim}${experiment}"
		cp "${pathPer}prof.inp.000" "${pathSim}${experiment}/prof.inp.001"
		cp "${pathPer}lscale.inp.000" "${pathSim}${experiment}/lscale.inp.001"
		cp "namoptions.${type}" "${pathSim}${experiment}/namoptions" 
	done
done
