#!/bin/bash
# Create folders for simulations and copy initial profiles
sigmax_array=( 000 002 004 008 016 )
sigmat_array=( 000 006 030 180 )
pathCases="../cases/"
pathOpenBC_synturb="${pathCases}openBC_synturb/"
mkdir -p $pathCases
mkdir -p $pathOpenBC_synturb
cd $pathOpenBC_synturb
for sigmat in ${sigmat_array[@]}; do
	for sigmax in ${sigmax_array[@]}; do
		experiment="x${sigmax}y${sigmax}z000t${sigmat}"
		# Skip simulation without smoothing as it's the same as the openBC simulation
    if [ $experiment = "x000y000z000t000" ]; then
      continue
    fi
		# create simulation directory
		mkdir -p $experiment
		cd $experiment
		rm -f dales4.4
		rm -f merge.sh
		rm -f openboundaries.inp.001.nc
		cp ../../../input/initial_profiles/prof.inp.xxx prof.inp.001
		cp ../../../input/initial_profiles/lscale.inp.xxx lscale.inp.001
		cp ../../../input/namoptions/namoptions.openBC_synturb namoptions
		ln -s "../../../input/boundary_input/openboundaries.inp.${experiment}.nc" openboundaries.inp.001.nc
		ln -s ../../../build/src/dales4.4 .
		ln -s ../../../scripts/merge.sh .
		cd ..
	done
done
cd ../../scripts/
