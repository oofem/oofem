nltrans_nonlin.out.tm
Quadrilateral element subjected to heat flux (Newton b.c)
NlTransientTransportProblem nsteps 3 deltat 3600 rtol 1.e-5 alpha 0.5 lumpedcapa exportfields 1 5 nmodules 0
#vtkxml tstep_all domain_all primvars 1 6
domain heattransfer
OutputManager tstep_all dofman_all element_all
ndofman 4 nelem 1 ncrosssect 1 nmat 1 nbc 1 nic 1 nltf 1
node 1 coords 3 0.00 0.00 0.00 ic 1 1
node 2 coords 3 0.04 0.00 0.00 ic 1 1
node 3 coords 3 0.00 0.12 0.00 ic 1 1
node 4 coords 3 0.04 0.12 0.00 ic 1 1
quad1ht 1 nodes 4 1 2 4 3 crosssect 1 mat 1 boundaryLoads 4 1 1 1 3
SimpleCS 1 thick 1.0
IsoHeat 1 d 2400. k 1.5 c 800.0
constantedgeload 1 loadTimeFunction 1 ndofs 1 components 1 -320.0 loadtype 2
InitialCondition 1 Conditions 1 u 0.0
ConstantFunction 1 f(t) 1.0
#%BEGIN_CHECK%
#NODE tStep 0.0 number 1 dof 1 unknown f value 10.0
#NODE tStep 0.0 number 3 dof 1 unknown f value 10.0
#NODE tStep 3600 number 1 dof 1 unknown f value 20.0
#NODE tStep 3600 number 3 dof 1 unknown f value 20.0
#NODE tStep 7200 number 1 dof 1 unknown f value 30.0
#NODE tStep 7200 number 3 dof 1 unknown f value 30.0
#%END_CHECK%

