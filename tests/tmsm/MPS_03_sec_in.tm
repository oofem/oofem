MPS_03_sec_out.tm
tests on creep at variable temperature and relative humidity
#
nltransienttransportproblem nsteps 26 alpha 0.5 rtol 1.e-10 lumpedcapa nsmax 1000 exportfields 2 5 6 prescribedtimes 26 1.e-10 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1. 2. 5. 10. 20. 50. 100. 200. 500. 1000. 2000. 5000. 10000.
# nmodules 1
#
# vtkxml tstep_all domain_all primvars 2 3 6
#
domain HeMa1
#
OutputManager tstep_all dofman_all element_all
ndofman 24 nelem 6 ncrosssect 1 nmat 1 nbc 5 nic 2 nltf 4
#
#
# NODES
#
node   1   coords 3  0.0  0.0  0.0 ic 2 1 2 bc 2 1 2
node   2   coords 3  0.1  0.0  0.0 ic 2 1 2 bc 2 1 2
node   3   coords 3  0.0  0.1  0.0 ic 2 1 2 bc 2 1 2
node   4   coords 3  0.1  0.1  0.0 ic 2 1 2 bc 2 1 2
#
node   5   coords 3  0.0  0.2  0.0 ic 2 1 2 bc 2 1 3
node   6   coords 3  0.1  0.2  0.0 ic 2 1 2 bc 2 1 3
node   7   coords 3  0.0  0.3  0.0 ic 2 1 2 bc 2 1 3
node   8   coords 3  0.1  0.3  0.0 ic 2 1 2 bc 2 1 3
#
node   9   coords 3  0.0  0.4  0.0 ic 2 1 2 bc 2 4 2
node   10  coords 3  0.1  0.4  0.0 ic 2 1 2 bc 2 4 2
node   11  coords 3  0.0  0.5  0.0 ic 2 1 2 bc 2 4 2
node   12  coords 3  0.1  0.5  0.0 ic 2 1 2 bc 2 4 2
#
node   13  coords 3  0.0  0.6  0.0 ic 2 1 2 bc 2 5 2
node   14  coords 3  0.1  0.6  0.0 ic 2 1 2 bc 2 5 2
node   15  coords 3  0.0  0.7  0.0 ic 2 1 2 bc 2 5 2
node   16  coords 3  0.1  0.7  0.0 ic 2 1 2 bc 2 5 2
#
node   17  coords 3  0.0  0.8  0.0 ic 2 1 2 bc 2 5 2
node   18  coords 3  0.1  0.8  0.0 ic 2 1 2 bc 2 5 2
node   19  coords 3  0.0  0.9  0.0 ic 2 1 2 bc 2 5 2
node   20  coords 3  0.1  0.9  0.0 ic 2 1 2 bc 2 5 2
#
node   21  coords 3  0.0  1.0  0.0 ic 2 1 2 bc 2 5 3
node   22  coords 3  0.1  1.0  0.0 ic 2 1 2 bc 2 5 3
node   23  coords 3  0.0  1.1  0.0 ic 2 1 2 bc 2 5 3
node   24  coords 3  0.1  1.1  0.0 ic 2 1 2 bc 2 5 3
#
#
# ELEMENTS
# reference basic + fixed temperature
quad1hmt   1   nodes 4   1 2 4 3 crossSect 1 mat 1 
#
# reference drying + fixed temperature
quad1hmt   2   nodes 4   5 6 8 7 crossSect 1 mat 1 
#
# elevated temperature + sealed
quad1hmt   3   nodes 4   9 10 12 11 crossSect 1 mat 1
#
# temperature cycle + sealed 
quad1hmt   4   nodes 4   13 14 16 15 crossSect 1 mat 1
#
# temperature cycle + sealed (temperature cycle damping)
quad1hmt   5   nodes 4   17 18 20 19 crossSect 1 mat 1
#
# temperature cycle + drying
quad1hmt   6   nodes 4   21 22 24 23 crossSect 1 mat 1 
#
# CROSSECTION
#
SimpleCS 1 thick 1.0 width 1.0
#
#
# MATERIAL
#
hemobaznajmat 1 d 2400. capa 1. C1 25.e-06 alpha0 0.05 hC 0.8 n 15. c 1.e6 k 1.
#
#
# BOUNDARY AND INITIAL CONDITIONS
#
# constant T - "room"
BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 298.
# constant h - "sealed"
BoundaryCondition 2 loadTimeFunction 2 prescribedvalue 1.
# drying
BoundaryCondition 3 loadTimeFunction 3 prescribedvalue 1.
# elevated temperature
BoundaryCondition 4 loadTimeFunction 1 prescribedvalue 348.
# temperature cycle
BoundaryCondition 5 loadTimeFunction 4 prescribedvalue 1.
#
InitialCondition 1 Conditions 1 u 298.
InitialCondition 2 Conditions 1 u 1.0
#
#
# TIME FUNCTION
#
ConstantFunction 1 f(t) 1.0
ConstantFunction 2 f(t) 1.0
#PiecewiseLinfunction 3 npoints  8 t 8 -1.e6 0.0  0.1 1. 10. 100. 1000. 10000. f(t) 8 1. 1. 1. 0.9 0.8 0.7 0.6 0.5
PiecewiseLinfunction 3 npoints  8 t 8 -8.64e10 0.0  8.64e3 8.64e4 8.64e5 8.64e6 8.64e7 8.64e8 f(t) 8 1. 1. 1. 0.9 0.8 0.7 0.6 0.5
#		       		      	    	     	    	       	      	     	      	            -1.e6  0.0  0.1   1.  10.  100. 200. 500. 1000. 2000. 5000. 10000.
#PiecewiseLinfunction 4 npoints  12 t 12 -1.e6 0.0  0.1 1. 10. 100. 200. 500. 1000. 2000. 5000. 10000. f(t) 12 298. 298. 298. 298. 300. 330. 300. 250. 250.  330.  300.  330.
PiecewiseLinfunction 4 npoints  12 t 12 -8.64000E+10 0.00000E+0 8.64000E+3 8.64000E+4 8.64000E+5 8.64000E+6 1.72800E+7 4.32000E+7 8.64000E+7 1.72800E+8 4.32000E+8 8.64000E+8  f(t) 12 298. 298. 298. 298. 300. 330. 300. 250. 250.  330.  300.  330.