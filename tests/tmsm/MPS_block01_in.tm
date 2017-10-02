MPS_block01_in.out.tm
OOFEM input file for thermal analysis, automatically generated, test for ConTemp, a brick with size x=0.6 y=0.5 z=0.5
TransientTransport nsteps 3 deltat 36000 rtolv 1.e-4 alpha 0.6 lumped exportfields 1 5 nmodules 1
#TransientTransport nsteps 3 deltat 36000 rtolf 1.e-6 alpha 0.6 lumped exportfields 1 5 nmodules 1
#vtkxml tstep_all domain_all primvars 1 6 vars 3 39 56 95 stype 2 timescale 2.77777e-4
errorcheck
domain HeatTransfer
OutputManager tstep_all dofman_output { 8 }
ndofman 27 nelem 8 ncrosssect 1 nmat 1 nbc 2 nic 1 nltf 7 nset 3
node 1 coords 3 0 0 0
node 2 coords 3 0.3 0 0
node 3 coords 3 0.6 0 0
node 4 coords 3 0 0.25 0
node 5 coords 3 0.3 0.25 0
node 6 coords 3 0.6 0.25 0
node 7 coords 3 0 0.5 0
node 8 coords 3 0.3 0.5 0
node 9 coords 3 0.6 0.5 0
node 10 coords 3 0 0 0.25
node 11 coords 3 0.3 0 0.25
node 12 coords 3 0.6 0 0.25
node 13 coords 3 0 0.25 0.25
node 14 coords 3 0.3 0.25 0.25
node 15 coords 3 0.6 0.25 0.25
node 16 coords 3 0 0.5 0.25
node 17 coords 3 0.3 0.5 0.25
node 18 coords 3 0.6 0.5 0.25
node 19 coords 3 0 0 0.5
node 20 coords 3 0.3 0 0.5
node 21 coords 3 0.6 0 0.5
node 22 coords 3 0 0.25 0.5
node 23 coords 3 0.3 0.25 0.5
node 24 coords 3 0.6 0.25 0.5
node 25 coords 3 0 0.5 0.5
node 26 coords 3 0.3 0.5 0.5
node 27 coords 3 0.6 0.5 0.5
Brick1ht 1 nodes 8  10 13 14 11 1 4 5 2
Brick1ht 2 nodes 8  11 14 15 12 2 5 6 3
Brick1ht 3 nodes 8  13 16 17 14 4 7 8 5
Brick1ht 4 nodes 8  14 17 18 15 5 8 9 6
Brick1ht 5 nodes 8  19 22 23 20 10 13 14 11
Brick1ht 6 nodes 8  20 23 24 21 11 14 15 12
Brick1ht 7 nodes 8  22 25 26 23 13 16 17 14
Brick1ht 8 nodes 8  23 26 27 24 14 17 18 15
Set 1 elementranges {(1 8)}
Set 2 elementboundaries 8  1 6  2 6  5 6  6 6
Set 3 elementboundaries 8  5 1  6 1  7 1  8 1
SimpleTransportCS 1 mat 1 set 1
HydratingConcreteMat 1 d 2445.52 k 1.52188 c 935.851 hydrationmodeltype 2 Qpot 471.15 masscement 199.956 b1 0.0001624 b2 0.0014 eta 7 dohinf 0.85 activationenergy 38300 castingTime 0. DoH1 0 P1 0
ConstantSurfaceLoad 1 loadTimeFunction 1 components 1 20 properties 1 a 19 loadtype 3 set 2
ConstantSurfaceLoad 2 loadTimeFunction 1 components 1 20 properties 1 a 19 loadtype 3 set 3
InitialCondition 1 conditions 1 u 20 dofs 1 10 set 1
ConstantFunction 1 f(t) 1.0
UsrDefLTF 2 f(t) 0+(0)*sin(2*3.14159*(t+28800-28800)/86400.)+h(-1.e+30*sin(2*3.14159*(t+28800-28800)/86400.))*(0)*sin(2*3.14159*(t+28800-28800)/86400.)
UsrDefLTF 3 f(t) 0+(0)*sin(2*3.14159*(t+28800-28800)/86400.)+h(-1.e+30*sin(2*3.14159*(t+28800-28800)/86400.))*(0)*sin(2*3.14159*(t+28800-28800)/86400.)
UsrDefLTF 4 f(t) 20+(0)*sin(2*3.14159*(t+28800-28800)/86400.)+h(-1.e+30*sin(2*3.14159*(t+28800-28800)/86400.))*(0)*sin(2*3.14159*(t+28800-28800)/86400.)
UsrDefLTF 5 f(t) 0+(0)*sin(2*3.14159*(t+28800-28800)/86400.)+h(-1.e+30*sin(2*3.14159*(t+28800-28800)/86400.))*(0)*sin(2*3.14159*(t+28800-28800)/86400.)
UsrDefLTF 6 f(t) 0+(0)*sin(2*3.14159*(t+28800-28800)/86400.)+h(-1.e+30*sin(2*3.14159*(t+28800-28800)/86400.))*(0)*sin(2*3.14159*(t+28800-28800)/86400.)
UsrDefLTF 7 f(t) 20+(0)*sin(2*3.14159*(t+28800-28800)/86400.)+h(-1.e+30*sin(2*3.14159*(t+28800-28800)/86400.))*(0)*sin(2*3.14159*(t+28800-28800)/86400.)
#
#%BEGIN_CHECK% tolerance 1e-6
#NODE tStep 10  number 8 dof 10 unknown d value 3.16472144e+01 tolerance 1e-4
#NODE tStep 15  number 8 dof 10 unknown d value 3.12574472e+01 tolerance 1e-4
#%END_CHECK%

