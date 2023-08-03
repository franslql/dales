#!/bin/bash
# Create folders for simulations and copy initial profiles
sigmax_array=( 000 002 004 008 016 )
sigmat_array=( 000 006 030 180 )
pathDALES="../build/src/"
pathInput="../input/"
pathCases="../cases/"
pathOpenBC="${pathCases}openBC/"
mkdir -p $pathCases
mkdir -p $pathOpenBC
for sigmat in ${sigmat_array[@]}; do
	for sigmax in ${sigmax_array[@]}; do
		experiment="x${sigmax}y${sigmax}z000t${sigmat}"
		# create simulation directory
		mkdir -p "${pathOpenBC}${experiment}"
		cp "${pathInput}initial_profiles/prof.inp.xxx" "${pathOpenBC}${experiment}/prof.inp.001"
		cp "${pathInput}initial_profiles/lscale.inp.xxx" "${pathOpenBC}${experiment}/lscale.inp.001"
		cp "${pathInput}namoptions/namoptions.openBC" "${pathOpenBC}${experiment}/namoptions"
		ln -s "${pathInput}boundary_input/openboundaries.inp.${experiment}.nc" "${pathOpenBC}${experiment}/openboundaries.inp.001.nc"
		ln -s "${pathDALES}dales4.4" "${pathOpenBC}${experiment}/."
		ln -s merge.sh "${pathOpenBC}${experiment}/." 
	done
done
