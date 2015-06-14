MPS_02_sec_out.tm
element 1 constant rh = 0.98, element 2 constant rh = 0.5, remaining elements drying from rh = 0.98 to 0.5
#
nltransienttransportproblem nsteps 26 alpha 0.5 rtol 1.e-10 lumpedcapa nsmax 1000 exportfields 1 6 prescribedtimes 26 1.e-10 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1. 2. 5. 10. 20. 50. 100. 200. 500. 1000. 2000. 5000. 10000.
# nmodules 1
#
# vtkxml tstep_all domain_all primvars 1 3
#
domain mass1transfer
#
OutputManager tstep_all dofman_all element_all
ndofman 24 nelem 6 ncrosssect 1 nmat 1 nbc 3 nic 2 nltf 3
#
#
# NODES
#
node   1   coords 3  0.0  0.0  0.0 ic 1 1 bc 1 1
node   2   coords 3  0.1  0.0  0.0 ic 1 1 bc 1 1
node   3   coords 3  0.0  0.1  0.0 ic 1 1 bc 1 1
node   4   coords 3  0.1  0.1  0.0 ic 1 1 bc 1 1
#
node   5   coords 3  0.0  0.2  0.0 ic 1 2 bc 1 2
node   6   coords 3  0.1  0.2  0.0 ic 1 2 bc 1 2
node   7   coords 3  0.0  0.3  0.0 ic 1 2 bc 1 2
node   8   coords 3  0.1  0.3  0.0 ic 1 2 bc 1 2
#
node   9   coords 3  0.0  0.4  0.0 ic 1 1 bc 1 3
node   10  coords 3  0.1  0.4  0.0 ic 1 1 bc 1 3
node   11  coords 3  0.0  0.5  0.0 ic 1 1 bc 1 3
node   12  coords 3  0.1  0.5  0.0 ic 1 1 bc 1 3
#
node   13  coords 3  0.0  0.6  0.0 ic 1 1 bc 1 3
node   14  coords 3  0.1  0.6  0.0 ic 1 1 bc 1 3
node   15  coords 3  0.0  0.7  0.0 ic 1 1 bc 1 3
node   16  coords 3  0.1  0.7  0.0 ic 1 1 bc 1 3
#
node   17  coords 3  0.0  0.8  0.0 ic 1 1 bc 1 3
node   18  coords 3  0.1  0.8  0.0 ic 1 1 bc 1 3
node   19  coords 3  0.0  0.9  0.0 ic 1 1 bc 1 3
node   20  coords 3  0.1  0.9  0.0 ic 1 1 bc 1 3
#
node   21  coords 3  0.0  1.0  0.0 ic 1 1 bc 1 3
node   22  coords 3  0.1  1.0  0.0 ic 1 1 bc 1 3
node   23  coords 3  0.0  1.1  0.0 ic 1 1 bc 1 3
node   24  coords 3  0.1  1.1  0.0 ic 1 1 bc 1 3
#
#
# ELEMENTS
#
quad1mt   1   nodes 4   1 2 4 3 crossSect 1 mat 1 
#
quad1mt   2   nodes 4   5 6 8 7 crossSect 1 mat 1 
#
quad1mt   3   nodes 4   9 10 12 11 crossSect 1 mat 1
#
quad1mt   4   nodes 4   13 14 16 15 crossSect 1 mat 1
#
quad1mt   5   nodes 4   17 18 20 19 crossSect 1 mat 1
#
quad1mt   6   nodes 4   21 22 24 23 crossSect 1 mat 1 
#
# CROSSECTION
#
SimpleCS 1 thick 1.0 width 1.0
#
#
# MATERIAL
#
isolinmoisturemat 1 d 2400. perm 1.e-5 capa 1.
#
#
# BOUNDARY & INITIAL CONDITIONS
#
BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.98
BoundaryCondition 2 loadTimeFunction 2 prescribedvalue 1.
BoundaryCondition 3 loadTimeFunction 3 prescribedvalue 1.
#
InitialCondition 1 Conditions 1 u 0.98
InitialCondition 2 Conditions 1 u 0.5
#
#
# TIME FUNCTION
#
ConstantFunction 1 f(t) 1.0
ConstantFunction 2 f(t) 0.5
PiecewiseLinfunction 3 npoints  8 t 8 -8.64e10 0.0  8.64e3 8.64e4 8.64e5 8.64e6 8.64e7 8.64e8 f(t) 8 0.98 0.98 0.98 0.9 0.8 0.7 0.6 0.5