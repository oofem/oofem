MPS_bench02_out.tm
OOFEM input file for thermal analysis, automatically generated, test for ConTemp, a brick with size x=3 y=0.5 z=0.5
TransientTransport nsteps 10 deltat 7200 rtolv 1.e-4 alpha 0.6 lstype 1 lstol 1e-5 lsiter 200 lsprecond 1 renumber 1 exportfields 1 5 nmodules 1
#TransientTransport nsteps 10 deltat 7200 alpha 0.6 exportfields 1 5 rtolf 1.e-5 lstype 1 lstol 1e-5 lsiter 200 lsprecond 1 nmodules 1
#vtkxml tstep_all domain_all primvars 1 6 vars 3 39 56 95 stype 2 timescale 2.77777e-4
errorcheck
domain HeatTransfer
OutputManager tstep_all dofman_output { 21 }
ndofman 72 nelem 28 ncrosssect 1 nmat 1 nbc 7 nic 1 nltf 7 nset 2
node 1 coords 3 0 0 0
node 2 coords 3 0.428571 0 0
node 3 coords 3 0.857143 0 0
node 4 coords 3 1.28571 0 0
node 5 coords 3 1.71429 0 0
node 6 coords 3 2.14286 0 0
node 7 coords 3 2.57143 0 0
node 8 coords 3 3 0 0
node 9 coords 3 0 0.25 0
node 10 coords 3 0.428571 0.25 0
node 11 coords 3 0.857143 0.25 0
node 12 coords 3 1.28571 0.25 0
node 13 coords 3 1.71429 0.25 0
node 14 coords 3 2.14286 0.25 0
node 15 coords 3 2.57143 0.25 0
node 16 coords 3 3 0.25 0
node 17 coords 3 0 0.5 0
node 18 coords 3 0.428571 0.5 0
node 19 coords 3 0.857143 0.5 0
node 20 coords 3 1.28571 0.5 0
node 21 coords 3 1.71429 0.5 0
node 22 coords 3 2.14286 0.5 0
node 23 coords 3 2.57143 0.5 0
node 24 coords 3 3 0.5 0
node 25 coords 3 0 0 0.25
node 26 coords 3 0.428571 0 0.25
node 27 coords 3 0.857143 0 0.25
node 28 coords 3 1.28571 0 0.25
node 29 coords 3 1.71429 0 0.25
node 30 coords 3 2.14286 0 0.25
node 31 coords 3 2.57143 0 0.25
node 32 coords 3 3 0 0.25
node 33 coords 3 0 0.25 0.25
node 34 coords 3 0.428571 0.25 0.25
node 35 coords 3 0.857143 0.25 0.25
node 36 coords 3 1.28571 0.25 0.25
node 37 coords 3 1.71429 0.25 0.25
node 38 coords 3 2.14286 0.25 0.25
node 39 coords 3 2.57143 0.25 0.25
node 40 coords 3 3 0.25 0.25
node 41 coords 3 0 0.5 0.25
node 42 coords 3 0.428571 0.5 0.25
node 43 coords 3 0.857143 0.5 0.25
node 44 coords 3 1.28571 0.5 0.25
node 45 coords 3 1.71429 0.5 0.25
node 46 coords 3 2.14286 0.5 0.25
node 47 coords 3 2.57143 0.5 0.25
node 48 coords 3 3 0.5 0.25
node 49 coords 3 0 0 0.5
node 50 coords 3 0.428571 0 0.5
node 51 coords 3 0.857143 0 0.5
node 52 coords 3 1.28571 0 0.5
node 53 coords 3 1.71429 0 0.5
node 54 coords 3 2.14286 0 0.5
node 55 coords 3 2.57143 0 0.5
node 56 coords 3 3 0 0.5
node 57 coords 3 0 0.25 0.5
node 58 coords 3 0.428571 0.25 0.5
node 59 coords 3 0.857143 0.25 0.5
node 60 coords 3 1.28571 0.25 0.5
node 61 coords 3 1.71429 0.25 0.5
node 62 coords 3 2.14286 0.25 0.5
node 63 coords 3 2.57143 0.25 0.5
node 64 coords 3 3 0.25 0.5
node 65 coords 3 0 0.5 0.5
node 66 coords 3 0.428571 0.5 0.5
node 67 coords 3 0.857143 0.5 0.5
node 68 coords 3 1.28571 0.5 0.5
node 69 coords 3 1.71429 0.5 0.5
node 70 coords 3 2.14286 0.5 0.5
node 71 coords 3 2.57143 0.5 0.5
node 72 coords 3 3 0.5 0.5
Brick1ht 1 nodes 8  25 33 34 26 1 9 10 2 BoundaryLoads 2 4 6
Brick1ht 2 nodes 8  26 34 35 27 2 10 11 3 BoundaryLoads 2 4 6
Brick1ht 3 nodes 8  27 35 36 28 3 11 12 4 BoundaryLoads 2 4 6
Brick1ht 4 nodes 8  28 36 37 29 4 12 13 5 BoundaryLoads 2 4 6
Brick1ht 5 nodes 8  29 37 38 30 5 13 14 6 BoundaryLoads 2 4 6
Brick1ht 6 nodes 8  30 38 39 31 6 14 15 7 BoundaryLoads 2 4 6
Brick1ht 7 nodes 8  31 39 40 32 7 15 16 8 BoundaryLoads 4 3 5 4 6
Brick1ht 8 nodes 8  33 41 42 34 9 17 18 10
Brick1ht 9 nodes 8  34 42 43 35 10 18 19 11
Brick1ht 10 nodes 8  35 43 44 36 11 19 20 12
Brick1ht 11 nodes 8  36 44 45 37 12 20 21 13
Brick1ht 12 nodes 8  37 45 46 38 13 21 22 14
Brick1ht 13 nodes 8  38 46 47 39 14 22 23 15
Brick1ht 14 nodes 8  39 47 48 40 15 23 24 16 BoundaryLoads 2 3 5
Brick1ht 15 nodes 8  49 57 58 50 25 33 34 26 BoundaryLoads 4 4 6 7 1
Brick1ht 16 nodes 8  50 58 59 51 26 34 35 27 BoundaryLoads 4 4 6 7 1
Brick1ht 17 nodes 8  51 59 60 52 27 35 36 28 BoundaryLoads 4 4 6 7 1
Brick1ht 18 nodes 8  52 60 61 53 28 36 37 29 BoundaryLoads 4 4 6 7 1
Brick1ht 19 nodes 8  53 61 62 54 29 37 38 30 BoundaryLoads 4 4 6 7 1
Brick1ht 20 nodes 8  54 62 63 55 30 38 39 31 BoundaryLoads 4 4 6 7 1
Brick1ht 21 nodes 8  55 63 64 56 31 39 40 32 BoundaryLoads 6 3 5 4 6 7 1
Brick1ht 22 nodes 8  57 65 66 58 33 41 42 34 BoundaryLoads 2 7 1
Brick1ht 23 nodes 8  58 66 67 59 34 42 43 35 BoundaryLoads 2 7 1
Brick1ht 24 nodes 8  59 67 68 60 35 43 44 36 BoundaryLoads 2 7 1
Brick1ht 25 nodes 8  60 68 69 61 36 44 45 37 BoundaryLoads 2 7 1
Brick1ht 26 nodes 8  61 69 70 62 37 45 46 38 BoundaryLoads 2 7 1
Brick1ht 27 nodes 8  62 70 71 63 38 46 47 39 BoundaryLoads 2 7 1
Brick1ht 28 nodes 8  63 71 72 64 39 47 48 40 BoundaryLoads 4 3 5 7 1
Set 1 elementranges {(1 28)}
Set 2 nodes 9 1 9 17 25 33 41 49 57 65
SimpleTransportCS 1 mat 1 set 1
HydratingConcreteMat 1 d 2403.24 k 1.471 c 961.929 hydrationmodeltype 2 Qpot 471.15 masscement 279.904 b1 0.0001624 b2 0.0014 eta 7 dohinf 0.85 activationenergy 38300 castingTime 0. DoH1 0 P1 0
BoundaryCondition 1 loadTimeFunction 1 dofs 1 10 values 1 20 set 2
ConstantSurfaceLoad 2 loadTimeFunction 1 components 1 20 properties 1 a 0 loadtype 3
ConstantSurfaceLoad 3 loadTimeFunction 1 components 1 20 properties 1 a 15 loadtype 3
ConstantSurfaceLoad 4 loadTimeFunction 1 components 1 20 properties 1 a 15 loadtype 3
ConstantSurfaceLoad 5 loadTimeFunction 1 components 1 0 properties 1 a 0 loadtype 3
ConstantSurfaceLoad 6 loadTimeFunction 1 components 1 0 properties 1 a 0 loadtype 3
ConstantSurfaceLoad 7 loadTimeFunction 1 components 1 20 properties 1 a 15 loadtype 3
InitialCondition 1 dofs 1 10 conditions 1 u 20 set 1
ConstantFunction 1 f(t) 1.0
UsrDefLTF 2 f(t) 20+(0)*sin(2*3.14159*(t+28800-28800)/86400.)+h(-1.e+30*sin(2*3.14159*(t+28800-28800)/86400.))*(0)*sin(2*3.14159*(t+28800-28800)/86400.)
UsrDefLTF 3 f(t) 20+(0)*sin(2*3.14159*(t+28800-28800)/86400.)+h(-1.e+30*sin(2*3.14159*(t+28800-28800)/86400.))*(0)*sin(2*3.14159*(t+28800-28800)/86400.)
UsrDefLTF 4 f(t) 20+(0)*sin(2*3.14159*(t+28800-28800)/86400.)+h(-1.e+30*sin(2*3.14159*(t+28800-28800)/86400.))*(0)*sin(2*3.14159*(t+28800-28800)/86400.)
UsrDefLTF 5 f(t) 0+(0)*sin(2*3.14159*(t+28800-28800)/86400.)+h(-1.e+30*sin(2*3.14159*(t+28800-28800)/86400.))*(0)*sin(2*3.14159*(t+28800-28800)/86400.)
UsrDefLTF 6 f(t) 0+(0)*sin(2*3.14159*(t+28800-28800)/86400.)+h(-1.e+30*sin(2*3.14159*(t+28800-28800)/86400.))*(0)*sin(2*3.14159*(t+28800-28800)/86400.)
UsrDefLTF 7 f(t) 20+(0)*sin(2*3.14159*(t+28800-28800)/86400.)+h(-1.e+30*sin(2*3.14159*(t+28800-28800)/86400.))*(0)*sin(2*3.14159*(t+28800-28800)/86400.)
#%BEGIN_CHECK% tolerance 1e-5
#NODE tStep 2  number 21 dof 10 unknown d value 2.03239872e+01
#NODE tStep 4  number 21 dof 10 unknown d value 2.16197503e+01
#NODE tStep 9  number 21 dof 10 unknown d value 3.25683074e+01
#NODE tStep 14  number 21 dof 10 unknown d value 3.95453679e+01
#NODE tStep 19  number 21 dof 10 unknown d value 3.99943586e+01
#%END_CHECK%
