# Piecewise Potential Vorticity Inversion (PPVI) tool

# Overview
This package ingests model forecast and analysis from the GFS, IFS and GFDL SHiELD models. The total potential vorticity (PV) field is calculated from total $\phi$, T, and u-,v-wind fields. The piecewise decomposition is performed which separates the total fields into the basic field (TC) and the perturbation field (environment). Then, the perturbation PV field is systematically divided into six regions emulating semi-permanent large-scale pressure systems: Bermuda High, Continental High and 4 quadrants in the ordinal directions (NE,SE,SW,NW). The PV in each of these regions is inverted to retrieve the individual mass $\phi$ and wind $\psi$ fields. The result is being able to quantify the steering flows from individual pressure systems that drive hurricane movement. 

# Acknowledgements
The original code is provided by Professor Chun-Chieh Wu's group at National Taiwan University and has been modified by Jan-Huey Chen (GFDL) and Tyler Barbero (CSU).
This package is written in Fortran 90, Python and C-shell wrapper scripts.

### Correspondence
tyler.barbero@colostate.edu\
jan-huey.chen@noaa.gov
