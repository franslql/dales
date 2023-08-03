#!/bin/bash
# Create folders for simulations and copy initial profiles
pathCases="../cases/"
pathPeriodic="${pathCases}periodic/"
mkdir -p $pathCases
mkdir -p $pathPeriodic
cd $pathPeriodic
rm -f dales4.4
rm -f merge.sh
cp ../../input/initial_profiles/prof.inp.xxx prof.inp.000
cp ../../input/initial_profiles/lscale.inp.xxx lscale.inp.000
cp ../../input/namoptions/namoptions.periodic namoptions
ln -s ../../build/src/dales4.4 .
ln -s ../../scripts/merge.sh .
cd ../../scripts/