MPS_04_out.tm
tests on creep with damage
#
nltransienttransportproblem nsteps 26 alpha 0.5 rtol 1.e-10 lumpedcapa nsmax 1000 exportfields 1 6 prescribedtimes 26 1.e-10 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1. 2. 5. 10. 20. 50. 100. 200. 500. 1000. 2000. 5000. 10000.
# nmodules 1 
#
# vtkxml tstep_all domain_all primvars 1 3
#
domain mass1transfer
#
OutputManager tstep_all dofman_all element_all
ndofman 12 nelem 3 ncrosssect 1 nmat 1 nbc 1 nic 1 nltf 1
#
#
# NODES
#
node   1   coords 3  0.0  0.0  0.0 ic 1 1 bc 1 1
node   2   coords 3  0.1  0.0  0.0 ic 1 1 bc 1 1
node   3   coords 3  0.0  0.1  0.0 ic 1 1 bc 1 1
node   4   coords 3  0.1  0.1  0.0 ic 1 1 bc 1 1
#
node   5   coords 3  0.0  0.2  0.0 ic 1 1 bc 1 1
node   6   coords 3  0.1  0.2  0.0 ic 1 1 bc 1 1
node   7   coords 3  0.0  0.3  0.0 ic 1 1 bc 1 1
node   8   coords 3  0.1  0.3  0.0 ic 1 1 bc 1 1
#
node   9   coords 3  0.0  0.4  0.0 ic 1 1 bc 1 1
node   10  coords 3  0.1  0.4  0.0 ic 1 1 bc 1 1
node   11  coords 3  0.0  0.5  0.0 ic 1 1 bc 1 1
node   12  coords 3  0.1  0.5  0.0 ic 1 1 bc 1 1
#
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
# BOUNDARY AND INITIAL CONDITIONS
#
BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 1.0
#
InitialCondition 1 Conditions 1 u 0.98
#
#
# TIME FUNCTION
#
PiecewiseLinfunction 1 npoints  8 t 8 -1.e6 0.0  0.1 1. 10. 100. 1000. 10000. f(t) 8 0.98 0.98 0.98 0.9 0.8 0.7 0.6 0.5