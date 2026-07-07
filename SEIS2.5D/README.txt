SEIS2.5D is a software package to calculate dynamic displacements and strains resulting from imposed earthquake sources in a spherical geometry. 
It solves the viscoelastic constitutive equation and linear stress-strain relation using the spectral element method (SEM) on a three-dimensional elastic
structure, with imposed source and domain boundary conditions. In the spectral element domain, it produces time series of dynamic displacements and strains at
sets of points on Earth’s surface and vertical profiles, as well as specified `receiver’ locations. 

The program package is intended as a tool in modeling forward and backward seismic wavefields in the case that the stress-strain relation is linear and the
elastic structure is two-dimensional. Such 2D structures could correspond to a contrast in material properties across, e.g., a long strike slip fault or a long
convergent plate boundary with a dipping downgoing plate.

The package has two dependencies which must be obtained before installation:

(1) Suitesparse is the linear algebra package needed to compile SEIS2PT5D. It may be downloaded from http://faculty.cse.tamu.edu/davis/suitesparse.html. This package is due to Timothy A Davis: http://faculty.cse.tamu.edu/davis/publications.html. The package should be placed in the subdirectory SEIS2.5D/Suitesparse
(2) The OpenBLAS library may be obtained from https://osdn.net/projects/sfnet_openblas/releases/. It should be placed in the subdirectory SEIS2.5D/OpenBLAS-0.3.13

To use this software, please find instructions in DOC/SEIS2.5D-tutorial.docx

License: This project is in the public domain.
Disclaimer: This software is preliminary or provisional and is subject to revision.


