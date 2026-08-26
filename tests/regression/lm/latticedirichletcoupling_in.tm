latticedirichletcoupling_out.tm
tm part of lattice Dirichlet hydro-mechanical boundary coupling
TransientTransport nsteps 5 deltat 1.0 alpha 1. nmodules 1
errorcheck
domain 2dMassLatticeTransport
OutputManager tstep_all dofman_all element_all
ndofman 5 nelem 4 ncrosssect 1 nmat 1 nbc 1 nic 0 nltf 1 nset 2
node 1 coords 2 1.0 0.0
node 2 coords 2 2.0 1.0
node 3 coords 2 1.0 2.0
node 4 coords 2 0.0 1.0
node 5 coords 2 1.0 1.0
latticemt2d 1 nodes 2 1 5 dim 2 mat 1 thick 1.0 width 1.0 gpCoords 2 1.0 0.5 crackwidth 0
latticemt2d 2 nodes 2 5 2 dim 2 mat 1 thick 1.0 width 1.0 gpCoords 2 1.5 1.0 crackwidth 0
latticemt2d 3 nodes 2 3 5 dim 2 mat 1 thick 1.0 width 1.0 gpCoords 2 1.0 1.5 crackwidth 0
latticemt2d 4 nodes 2 4 5 dim 2 mat 1 thick 1.0 width 1.0 gpCoords 2 0.5 1.0 crackwidth 0
simpletransportcs 1 mat 1 set 1
latticetransmat 1 d 1.e-9 k 1.e-13 vis 1.e-9 thetas 0.0924 thetar 0. ctor 0.001
LatticeDirichletCoupling 1 loadTimeFunction 1 dofs 1 11 values 1 0. couplingelements 2 1 3 set 2
ConstantFunction 1 f(t) 1.
Set 1 elements 4 1 2 3 4
Set 2 nodes 1 5
#%BEGIN_CHECK% tolerance 1.e-5
#NODE tStep 2 number 4 dof 11 unknown d value -7.029050e+02
#NODE tStep 3 number 4 dof 11 unknown d value -1.405810e+03
#%END_CHECK%
