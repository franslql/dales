#!/bin/bash
# Create folders for simulations and copy input files
tauh_array=( 000 060 )
pbc_array=( 2 4 )
dint_array=( 00 05 )
pathCases="../cases/"
pathSensitivity="${pathCases}sensitivity/"
mkdir -p $pathCases
mkdir -p $pathSensitivity
cd $pathSensitivity
# Create simulations for variation in tauh
for tauh in ${tauh_array[@]}; do
  experiment="tauh${tauh}"
  mkdir -p $experiment
  # create simulation directory
  mkdir -p $experiment
  cd $experiment
  rm -f dales4.4
  rm -f merge.sh
  rm -f openboundaries.inp.001.nc
  cp ../../../input/initial_profiles/prof.inp.xxx prof.inp.001
  cp ../../../input/initial_profiles/lscale.inp.xxx lscale.inp.001
  cp "../../../input/namoptions/namoptions.openBC_tauh${tauh}" namoptions
  ln -s ../../../input/boundary_input/openboundaries.inp.x000y000z000t000.nc openboundaries.inp.001.nc
  ln -s ../../../build/src/dales4.4 .
  ln -s ../../../scripts/merge.sh .
  cd ..
done
# Create simulations for variation in pbc
for pbc in ${pbc_array[@]}; do
  experiment="pbc${pbc}"
  mkdir -p $experiment
  # create simulation directory
  mkdir -p $experiment
  cd $experiment
  rm -f dales4.4
  rm -f merge.sh
  rm -f openboundaries.inp.001.nc
  cp ../../../input/initial_profiles/prof.inp.xxx prof.inp.001
  cp ../../../input/initial_profiles/lscale.inp.xxx lscale.inp.001
  cp "../../../input/namoptions/namoptions.openBC_pbc${pbc}" namoptions
  ln -s ../../../input/boundary_input/openboundaries.inp.x000y000z000t000.nc openboundaries.inp.001.nc
  ln -s ../../../build/src/dales4.4 .
  ln -s ../../../scripts/merge.sh .
  cd ..
done
# Create simulations for variation in dint
for dint in ${dint_array[@]}; do
  experiment="dint${dint}"
  mkdir -p $experiment
  # create simulation directory
  mkdir -p $experiment
  cd $experiment
  rm -f dales4.4
  rm -f merge.sh
  rm -f openboundaries.inp.001.nc
  cp ../../../input/initial_profiles/prof.inp.xxx prof.inp.001
  cp ../../../input/initial_profiles/lscale.inp.xxx lscale.inp.001
  cp "../../../input/namoptions/namoptions.openBC_dint${dint}" namoptions
  ln -s ../../../input/boundary_input/openboundaries.inp.x000y000z000t000.nc openboundaries.inp.001.nc
  ln -s ../../../build/src/dales4.4 .
  ln -s ../../../scripts/merge.sh .
  cd ..
done
# Create simulation without buoyancy term at the top boundary
experiment="lbuoyoff"
mkdir -p $experiment
# create simulation directory
mkdir -p $experiment
cd $experiment
rm -f dales4.4
rm -f merge.sh
rm -f openboundaries.inp.001.nc
cp ../../../input/initial_profiles/prof.inp.xxx prof.inp.001
cp ../../../input/initial_profiles/lscale.inp.xxx lscale.inp.001
cp ../../../input/namoptions/namoptions.openBC_lbuoyoff namoptions
ln -s ../../../input/boundary_input/openboundaries.inp.x000y000z000t000.nc openboundaries.inp.001.nc
ln -s ../../../build/src/dales4.4 .
ln -s ../../../scripts/merge.sh .
cd ..
cd ../../scripts/