This software was originally developed by Duo Li (Li & Liu, 2016). This is an updated version developed for the Costa Rica SSE and earthquake sequence model.

Li, D., & Liu, Y. (2016). Spatiotemporal evolution of slow slip events in a nonplanar fault model for northern Cascadia subduction zone. Journal of Geophysical Research: Solid Earth, 121(9), 6828-6845.

How to use this package:
1. Compute the Green's function. Please refer to readme.txt in the ./GreenFunction directory.
2. Create the output directory:   mkdir h15_out
3. Compile the program:

   mpifort phy3d_module_v3-2.f90 3dtriSSE_v3-2_p3.f90 -o program_run

   This program is compiled with the Intel Fortran compiler.

This program uses 148 CPUs. If you want to change the number of CPUs, make sure the number of cells is divisible by the number of CPUs. Change the variable nproc in ./GreenFunction (check readme.txt in the directory). Edit the parameter.txt file:

line 9: 5 132904 898 148 100 8000 34225 6947 3190 77086 4792 !Nab,Nt_all,Nt,nprocs,hnucl,nre,nsouth,nnorth,neqzone,nssezone,nnshallow

Modify the third and fourth input, the third input is n_cell/nproc, the fourth input is nproc


