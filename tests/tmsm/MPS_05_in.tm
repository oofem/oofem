MPS_05_out.tm
tests on creep at variable temperature 
#
nltransienttransportproblem nsteps 26 alpha 0.5 rtol 1.e-10 lumpedcapa nsmax 1000 exportfields 1 5 prescribedtimes 26 1.e-10 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1. 2. 5. 10. 20. 50. 100. 200. 500. 1000. 2000. 5000. 10000.
# nmodules 1
#
# vtkxml tstep_all domain_all primvars 2 3 6
#
domain HeatTransfer
#
OutputManager tstep_all dofman_all element_all
ndofman 24 nelem 6 ncrosssect 1 nmat 1 nbc 4 nic 2 nltf 4
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
node   13  coords 3  0.0  0.6  0.0 ic 1 2 bc 1 4
node   14  coords 3  0.1  0.6  0.0 ic 1 2 bc 1 4
node   15  coords 3  0.0  0.7  0.0 ic 1 2 bc 1 4
node   16  coords 3  0.1  0.7  0.0 ic 1 2 bc 1 4
#
node   17  coords 3  0.0  0.8  0.0 ic 1 1 bc 1 1
node   18  coords 3  0.1  0.8  0.0 ic 1 1 bc 1 1
node   19  coords 3  0.0  0.9  0.0 ic 1 1 bc 1 1
node   20  coords 3  0.1  0.9  0.0 ic 1 1 bc 1 1
#
node   21  coords 3  0.0  1.0  0.0 ic 1 2 bc 1 2
node   22  coords 3  0.1  1.0  0.0 ic 1 2 bc 1 2
node   23  coords 3  0.0  1.1  0.0 ic 1 2 bc 1 2
node   24  coords 3  0.1  1.1  0.0 ic 1 2 bc 1 2
#
#
# ELEMENTS
# element 1: expansion - K
quad1ht   1   nodes 4   1 2 4 3 crossSect 1 mat 1 
#
# element 2: expansion - C
quad1ht   2   nodes 4   5 6 8 7 crossSect 1 mat 1 
#
# element 3: accelerated creep - K
quad1ht   3   nodes 4   9 10 12 11 crossSect 1 mat 1
#
# element 4: accelerated creep - C
quad1ht   4   nodes 4   13 14 16 15 crossSect 1 mat 1
#
# element 5: TTC - K
quad1ht   5   nodes 4   17 18 20 19 crossSect 1 mat 1
#
# element 6: TTC - C
quad1ht   6   nodes 4   21 22 24 23 crossSect 1 mat 1 
#
# CROSSECTION
#
SimpleCS 1 thick 1.0 width 1.0
#
#
# MATERIAL
#
isoheat 1 d 2400. c 1.e6 k 1.
#
#
# BOUNDARY AND INITIAL CONDITIONS
#
# temperature cycle - K
BoundaryCondition 1 loadTimeFunction 3 prescribedvalue 1.
# temperature cycle - C
BoundaryCondition 2 loadTimeFunction 4 prescribedvalue 1.
# elevated temperature - K
BoundaryCondition 3 loadTimeFunction 1 prescribedvalue 348.15
# elevated temperature - C
BoundaryCondition 4 loadTimeFunction 2 prescribedvalue 75.
#
InitialCondition 1 Conditions 1 u 298.15
InitialCondition 2 Conditions 1 u 25.
#
#
# TIME FUNCTION
#
ConstantFunction 1 f(t) 1.0
ConstantFunction 2 f(t) 1.0
#		       		      	    	     	    	       	      	     	      	            -1.e6  0.0  0.1   1.  10.  100. 200. 500. 1000. 2000. 5000. 10000.
PiecewiseLinfunction 3 npoints  12 t 12 -1.e6 0.0  0.1 1. 10. 100. 200. 500. 1000. 2000. 5000. 10000. f(t) 12 298.15 298.15 298.15 298.15 300.15 330.15 300.15 250.15 250.15  330.15  298.15  298.15
#		       		      	    	     	    	       	      	     	      	            -1.e6  0.0  0.1   1.  10.  100. 200. 500. 1000. 2000. 5000. 10000.
PiecewiseLinfunction 4 npoints  12 t 12 -1.e6 0.0  0.1 1. 10. 100. 200. 500. 1000. 2000. 5000. 10000. f(t) 12  25. 25.  25.   25. 27.   57.  27. -23. -23.   57.  25.   25.
#
#%BEGIN_CHECK%
#TIME
#NODE number 4 dof 10 unknown d
#NODE number 8 dof 10 unknown d
#NODE number 12 dof 10 unknown d
#NODE number 16 dof 10 unknown d
#NODE number 20 dof 10 unknown d
#NODE number 24 dof 10 unknown d
#%END_CHECK%