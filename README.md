# DALES - Dutch Atmospheric Large Eddy Simulation

[![DOI](https://zenodo.org/badge/32735454.svg)](https://zenodo.org/badge/latestdoi/32735454)

## Documentation
The following documents are included in the DALES repository (but are not fully up to date):

* [User manual](https://github.com/dalesteam/dales/blob/master/utils/doc/input/dales-manual.pdf)
* Documentation of the configuration options [Namoptions.pdf](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf)
* [readme.pdf](https://github.com/dalesteam/dales/blob/master/utils/doc/input/readme.pdf)
* [INSTALL.md](https://github.com/dalesteam/dales/blob/master/INSTALL.md) contains installation instructions
* [CHANGELOG.md](https://github.com/dalesteam/dales/blob/master/CHANGELOG.md) documents the changes between versions
* The [DALES Wiki](https://github.com/dalesteam/dales/wiki/) contains [Installation instructions](https://github.com/dalesteam/dales/wiki/Installation-notes) for various systems and notes for [visualizing and processing](https://github.com/dalesteam/dales/wiki/Visualizing-and-processing-output) the model output

## License
DALES is made available under the terms of the GNU GPL version 3, see the file COPYING for details.

## References
Formulation of the Dutch Atmospheric Large-Eddy Simulation (DALES) and overview of its applications, T. Heus et al, [Geosci. Model Dev., 3, 415-444, 2010](https://doi.org/10.5194/gmd-3-415-2010)

Ouwersloot, H.G., Moene, A.F., Attema, J.J. et al. Large-Eddy Simulation Comparison of Neutral Flow Over a Canopy: Sensitivities to Physical and Numerical Conditions, and Similarity to Other Representations. [Boundary-Layer Meteorol 162, 71–89 (2017)](https://doi.org/10.1007/s10546-016-0182-5)

## Recreating manuscript simulations and figures
Notes: need cdo and fftw3 libraries. Results will be in results/.
1. Build DALES (Done)
export SYST=gnu-fast
mkdir build
cd build/
cmake .. -DUSE_FFTW=True
make

2. Create file structure for periodic simulation
cd scripts/
./createSimulationsPeriodic.sh

3. Run periodic simulation
cd scripts/
./run_periodic.sh

4. Create file structure for openBC
cd scripts/
./createSimulationsOpenBC.sh

5. Create boundary input for openBC
cd scripts/
python createInputOpenBC.py

6. Run openBC simulations
cd scripts/
./run_openBC.sh

7. Visualise results openBC
cd scripts/
python visualizePerOpenBC.py
python visualizeCrossxy.py
python visualizeTkexz.py
python visualizeWavelet.py

8. Create file structure for openBC_synturb
cd scripts/
./createSimulationsOpenBC_synturb.sh

9. Create boundary input for openBC_synturb
cd scripts/
python createInputOpenBC_synturb.py

10. Run openBC_synturb simulations
cd scripts/
./runOpenBC_synturb.sh

11. Visualise results openBC_synturb
cd scripts/
python visualizeCrossxy_synturb.py
python visualizeTkexz_synturb.py
python visualizeWavelet_synturb.py

12. Create sensitivity simulations
cd scripts/
./createSimulationsSensitivity.sh

13. Run sensitivity simulations
cd scripts
./runSensitivity.sh

14. Visualise results sensitivity simulations
cd scripts/
python visualizeSensitivity.py