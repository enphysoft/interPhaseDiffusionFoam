# interPhaseDiffusionFoam
This solver, interPhaseDiffusionFoam, is to solve solute diffusion phenomena in two-phase flow, especially air and water. Passive convection of solute transport is implemented using the basic algorithm of interFoam under the scope of the volume of the fluid (VOF) algorithm, and this new algorithm simulates solute diffusion more accurately at the air-water interface.

## Installation
After downloading this package, go to the directory and do the followings:
1. log in as a sudo (in Linux only): $ sudo bash
2. set up an OpenFOAM environment (for OF5): $ source /opt/openfoam5/etc/bashrc
3. make the solver: $ wmake

Solely developed by Albert S. Kim
