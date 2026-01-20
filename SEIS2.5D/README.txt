SEIS2.5D is a software package to calculate dynamic displacements and strains resulting from imposed earthquake sources in a spherical geometry. It solves the viscoelastic constitutive equation and linear stress-strain relation using the spectral element method (SEM) on a three-dimensional elastic structure, with imposed source and domain boundary conditions. In the spectral element domain, it produces time series of dynamic displacements and strains at sets of points on Earth’s surface and vertical profiles, as well as specified `receiver’ locations. 
The program package is intended as a tool in modeling coseismic or postseismic deformation from seismic sources.

Folders OpenBLAS-0.3.13 and Suitesparse comprise the linear algebra package needed to compile VISCO3D.  Placeholders for these folders are in the VISCO3D folder.  OpenBLAS may be obtained from https://osdn.net/projects/sfnet_openblas/releases/
Suitesparse may be obtained as the package SuiteSparse v4.4.4 downloaded from 
http://faculty.cse.tamu.edu/davis/suitesparse.html
This package is due to Timothy A Davis:  http://faculty.cse.tamu.edu/davis/publications.html


To compile SEIS2.5D: 

1) Un-tar SEIS2.5D.tar which will create the VISCO3D folder

2) Obtain the OpenBLAS-0.3.13 and Suitesparse folders from the online sources described above

3) Go into the SEIS2.5D/MAINPROG folder 
    > make libs 
    > make umf4_f77zwrapper64.o

4) Edit the file SEIS2.5D/SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk

In the line that reads 
BLAS = /Users/fpollitz-local/fred4/CIDER2015/VISCO3D/OpenBLAS-0.3.13/libopenblas.a
Replace the right hand side with the path to that library on your machine.

4) Go into the SEIS2.5D/MAINPROG folder
    > make seis2pt5d

