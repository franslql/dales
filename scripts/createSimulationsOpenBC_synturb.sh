#!/bin/bash
# Create folders for simulations and copy initial profiles
sigmax_array=( 000 002 004 008 016 )
sigmat_array=( 000 006 030 180 )
pathDALES="../build/src/"
pathInput="../input/"
pathCases="../cases/"
pathOpenBC_synturb="${pathCases}openBC_synturb/"
mkdir -p $pathCases
mkdir -p $pathOpenBC_synturb
for sigmat in ${sigmat_array[@]}; do
	for sigmax in ${sigmax_array[@]}; do
		experiment="x${sigmax}y${sigmax}z000t${sigmat}"
		# Skip simulation without smoothing as it's the same as the openBC simulation
    if [ $experiment = "x000y000z000t000" ]; then
      continue
    fi
		# create simulation directory
		mkdir -p "${pathOpenBC_synturb}${experiment}"
		cp "${pathInput}profiles/prof.inp.xxx" "${pathOpenBC_synturb}${experiment}/prof.inp.001"
		cp "${pathInput}profiles/lscale.inp.xxx" "${pathOpenBC_synturb}${experiment}/lscale.inp.001"
		cp "${pathInput}namoptions/namoptions.openBC_synturb" "${pathOpenBC_synturb}${experiment}/namoptions"
		ln -s "${pathInput}boundary_input/openboundaries.inp.${experiment}.nc" "${pathOpenBC_synturb}${experiment}/openboundaries.inp.001.nc"
		ln -s "${pathDALES}dales4.4" "${pathOpenBC_synturb}${experiment}/."
		ln -s merge.sh "${pathOpenBC_synturb}${experiment}/." 
	done
done
