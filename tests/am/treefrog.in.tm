treefrog.tm.out
Additive manufacturing - heat transfer
# Note: for production run increase nsteps to 8600
TransientTransport nsteps 20 deltat 1.0 alpha 0.5 lumped exportfields 1 5 nmodules 1 smtype 4 lstype 1 lsprecond 4 lstol 1e-9
vtkxml tstep_step 10 domain_all primvars 1 6
domain heattransfer
OutputManager
ndofman 0 nelem 0 ncrosssect 1 nmat 3 nbc 0 nic 0 nltf 1 nset 0
# Do not set any sets here - done by voxel activator
SimpleTransportCS 1 thickness 1.0 mat 1
TwoPhaseMat 1 mat 2 2 3
# air
IsoHeat 2 d 1250. k 0.2 c 1300.0 
# plastic
IsoHeat 3 d 1250. k $0.0002*te+0.1708$ c $5.09*te+974$ 
ConstantFunction 1 f(t) 1.0
