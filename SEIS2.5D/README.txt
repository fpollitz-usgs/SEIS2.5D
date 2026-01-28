SEIS2.5D is a software package to calculate dynamic displacements and strains resulting from imposed earthquake sources in a spherical geometry. It solves the viscoelastic constitutive equation and linear stress-strain relation using the spectral element method (SEM) on a three-dimensional elastic structure, with imposed source and domain boundary conditions. In the spectral element domain, it produces time series of dynamic displacements and strains at sets of points on Earth’s surface and vertical profiles, as well as specified `receiver’ locations. 

The program package is intended as a tool in modeling forward and backward seismic wavefields in the case that the stress-strain relation is linear and the elastic structure is two-dimensional. Such 2D structures could correspond to a contrast in material properties across, e.g., a long strike slip fault or a long convergent plate boundary with a dipping downgoing plate.
The package has two dependencies which must be obtained before installation:
(1) Suitesparse is the linear algebra package needed to compile SEIS2PT5D. It may be downloaded from http://faculty.cse.tamu.edu/davis/suitesparse.html. This package is due to Timothy A Davis: http://faculty.cse.tamu.edu/davis/publications.html. The package should be placed in the subdirectory SEIS2.5D/Suitesparse
(2) The OpenBLAS library may be obtained from https://osdn.net/projects/sfnet_openblas/releases/. It should be placed in the subdirectory SEIS2.5D/OpenBLAS-0.3.13
To compile SEIS2PT5D:
	1	Un-tar SEIS2.5D.tar which will create the SEIS2.5D folder
	2	Go into the SEIS2.5D/MAINPROG folder 
        make libs 
        make umf4_f77zwrapper64.o 
	3	Edit the file SEIS2.5D/SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk
In the line that reads BLAS = [] Replace the right hand side with the path to that library on your machine.
	4	Go into the SEIS2.5D/MAINPROG folder make seis2pt5d 
To run the example computation:
	1	Go into the SEIS2.5D/EXAMPLES folder
	2	Have a working copy of seis2pt5d in this folder (e.g., copied from the MAINPROG folder)
	3	Run example1.x-2proc
The goal of the example is to calculate strain seismograms at all of the GLL nodes on a particular vertical plane plus three-component velocity seismograms at a set of receivers (the 10 receivers listed in receivers-latlondepOLYM.txt). This is done in three steps.
First, seis2pt5d is run with a first argument of 1, which instructs it to compute bookkeeping indices that pertain to the sparse matrices used in the spectral element model.
Second, seis2pt5d is run with a first argument of 0, which instructs it to compute the synthetic seismograms (various strain and velocity seismograms) frequency-domain, and it stores the frequency domain responses on hard disk in the Work1 and Work2 folders. The task of looping over frequency is divided equally between two CPUs -- the first handles the fraction 0 - 0.5 of the frequencies, the second handles the fraction 0.5 - 1.0 of the frequencies, hence the arguments in the standard input to seis2pt5d.
Third, seis2pt5d is run with a first argument of 4, which instructs it to gather the frequency-domain response files (rec-j-m1 displ-j-m1; rec-j-m2 displ-j-m2; vertp-j files, which were stored in the second step in the Work1 and Work2 folders) compute inverse FT's at a large set of GLL nodes on a vertical plane and at the 10 receiver points.
The output files seis2pt5d.outstrains_vertp and seis2pt5d.outxyz_rec are then renamed to seis2pt5d.outstrains_vertp-example1 and seis2pt5d.outxyz_rec-DR-example1, respectively).
License: This project is in the public domain.
Disclaimer: This software is preliminary or provisional and is subject to revision.

