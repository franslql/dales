#!/bin/bash
# Create folders for simulations and copy initial profiles
pathInput="../input/"
pathCases="../cases/"
pathPeriodic="${pathCases}periodic/"
mkdir $pathCases
mkdir $pathPeriodic
cp "${pathInput}profiles/prof.inp.xxx" "${pathPeriodic}${experiment}/prof.inp.000"
cp "${pathInput}profiles/lscale.inp.xxx" "${pathPeriodic}${experiment}/lscale.inp.000"
cp "${pathInput}namoptions/namoptions.periodic" "${pathPeriodic}${experiment}/namoptions"