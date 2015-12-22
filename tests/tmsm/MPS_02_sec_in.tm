MPS_02_sec.out.tm
element 1 constant rh = 0.98, element 2 constant rh = 0.5, remaining elements drying from rh = 0.98 to 0.5
#
nltransienttransportproblem nsteps 26 alpha 0.5 rtol 1.e-10 lumpedcapa nsmax 1000 exportfields 1 6 prescribedtimes 26 1.e-10 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1. 2. 5. 10. 20. 50. 100. 200. 500. 1000. 2000. 5000. 10000.
#transienttransport nsteps 26 alpha 0.5 rtolf 1.e-10 lumped maxiter 1000 exportfields 1 6 prescribedtimes 26 1.e-10 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1. 2. 5. 10. 20. 50. 100. 200. 500. 1000. 2000. 5000. 10000.
# nmodules 1
#
# vtkxml tstep_all domain_all primvars 1 3
#
domain mass1transfer
#
OutputManager tstep_all dofman_all element_all
ndofman 24 nelem 6 ncrosssect 1 nmat 1 nbc 3 nic 2 nltf 3 nset 5
#
#
# NODES
#
node   1   coords 3  0.0  0.0  0.0
node   2   coords 3  0.1  0.0  0.0
node   3   coords 3  0.0  0.1  0.0
node   4   coords 3  0.1  0.1  0.0
#
node   5   coords 3  0.0  0.2  0.0
node   6   coords 3  0.1  0.2  0.0
node   7   coords 3  0.0  0.3  0.0
node   8   coords 3  0.1  0.3  0.0
#
node   9   coords 3  0.0  0.4  0.0
node   10  coords 3  0.1  0.4  0.0
node   11  coords 3  0.0  0.5  0.0
node   12  coords 3  0.1  0.5  0.0
#
node   13  coords 3  0.0  0.6  0.0
node   14  coords 3  0.1  0.6  0.0
node   15  coords 3  0.0  0.7  0.0
node   16  coords 3  0.1  0.7  0.0
#
node   17  coords 3  0.0  0.8  0.0
node   18  coords 3  0.1  0.8  0.0
node   19  coords 3  0.0  0.9  0.0
node   20  coords 3  0.1  0.9  0.0
#
node   21  coords 3  0.0  1.0  0.0
node   22  coords 3  0.1  1.0  0.0
node   23  coords 3  0.0  1.1  0.0
node   24  coords 3  0.1  1.1  0.0
#
#
# ELEMENTS
#
quad1mt   1   nodes 4   1 2 4 3
#
quad1mt   2   nodes 4   5 6 8 7
#
quad1mt   3   nodes 4   9 10 12 11
#
quad1mt   4   nodes 4   13 14 16 15
#
quad1mt   5   nodes 4   17 18 20 19
#
quad1mt   6   nodes 4   21 22 24 23
#
Set 1 elementranges {(1 6)}
Set 2 elements 1 1
Set 3 elements 1 2
Set 4 elements 4 3 4 5 6
Set 5 elements 5 1 3 4 5 6
#
# CROSSECTION
#
SimpleTransportCS 1 thickness 1.0 mat 1 set 1
#
#
# MATERIAL
#
isolinmoisturemat 1 d 2400. perm 1.e-5 capa 1.
#
#
# BOUNDARY & INITIAL CONDITIONS
#
BoundaryCondition 1 loadTimeFunction 1 dofs 1 14 values 1 0.98 set 2
BoundaryCondition 2 loadTimeFunction 2 dofs 1 14 values 1 1. set 3
BoundaryCondition 3 loadTimeFunction 3 dofs 1 14 values 1 1. set 4
#
InitialCondition 1 Conditions 1 u 0.98 dofs 1 14 set 5
InitialCondition 2 Conditions 1 u 0.5 dofs 1 14 set 3
#
#
# TIME FUNCTION
#
ConstantFunction 1 f(t) 1.0
ConstantFunction 2 f(t) 0.5
PiecewiseLinfunction 3 npoints  8 t 8 -8.64e10 0.0  8.64e3 8.64e4 8.64e5 8.64e6 8.64e7 8.64e8 f(t) 8 0.98 0.98 0.98 0.9 0.8 0.7 0.6 0.5