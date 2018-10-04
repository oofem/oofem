concretefcmviscoelasticcasting_out.tm
tests on creep with damage
#
transienttransport nsteps 26 alpha 0.5 exportfields 1 6 prescribedtimes 26 1.e-10 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1. 2. 5. 10. 20. 50. 100. 200. 500. 1000. 2000. 5000. 10000.
# nmodules 1 
#
# vtkxml tstep_all domain_all primvars 1 3
#
domain mass1transfer
#
OutputManager tstep_all dofman_all element_all
ndofman 4 nelem 1 ncrosssect 1 nmat 1 nbc 1 nic 1 nltf 1 nset 1
#
#
# NODES
#
node   1   coords 3  0.0  0.0  0.0
node   2   coords 3  1.0  0.0  0.0
node   3   coords 3  0.0  1.0  0.0
node   4   coords 3  1.0  1.0  0.0
#
#
#
#
#
# ELEMENTS
#
quad1mt   1   nodes 4   1 2 4 3
#
#
#
Set 1 elements 1 1
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
# BOUNDARY AND INITIAL CONDITIONS
#
BoundaryCondition 1 loadTimeFunction 1 dofs 1 14 values 1 1.0 set 1
#
InitialCondition 1 conditions 1 u 0.98 dofs 1 14 set 1
#
#
# TIME FUNCTION
#
PiecewiseLinfunction 1 npoints  8 t 8 -1.e6 14.0  14.1 15. 24. 114. 1014. 10014. f(t) 8 0.98 0.98 0.98 0.9 0.8 0.7 0.6 0.5
