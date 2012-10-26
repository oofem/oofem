############################################################
#
# Beam structure example from OOFEM documentation
# built from python
#
############################################################
import oofemlib

# engngModel
problem = oofemlib.linearStatic(nSteps=3, outFile="test2.out")

# domain (if no engngModel specified to domain, it is asigned to the last one created)
domain = oofemlib.domain(1, 1, problem, oofemlib.domainType._2dBeamMode, tstep_all=True, dofman_all=True, element_all=True)
problem.setDomain(1, domain)
vtkxmlModule = oofemlib.vtkxml(1,problem,tstep_all=True,vars=[1,4],primvars=[1])


# time functions (time functions parameter f(t) has to be written as f_t - due to the parenthesis, otherwise all kw keys match oofem input file ones)
# if no domain specified, ltf is asigned to the last one created (valid for all FEMComponents classes)
# if no number is specified, domain.umberOfLoadTimeFunctions+1 is used (valid for all FEMComponents classes, e.g. domain.numberOfElements+1 is used for elements instead ot numberOfLoadTimeFunctions)
ltf1 = oofemlib.peakFunction(1, domain, t=1, f_t=1)
ltf2 = oofemlib.peakFunction(2, domain, t=2, f_t=1)
ltf3 = oofemlib.peakFunction(3, domain, t=3, f_t=1)
ltfs = (ltf1, ltf2, ltf3)

# boundary conditions
# loadTimeFunction parameter can be specified as int value or as LoadTimeFunction itself (valid for all objects with giveNumber() method)
bc1   = oofemlib.boundaryCondition(    1, domain, loadTimeFunction=1,    prescribedValue=0.0)
bc2   = oofemlib.boundaryCondition(    2, domain, loadTimeFunction=2,    prescribedValue=-.006e-3)
eLoad = oofemlib.constantEdgeLoad(     3, domain, loadTimeFunction=ltf1, components=(0.,10.,0.), loadType=3, ndofs=3)
nLoad = oofemlib.nodalLoad(            4, domain, loadTimeFunction=1,    components=(-18.,24.,0.))
tLoad = oofemlib.structTemperatureLoad(5, domain, loadTimeFunction=ltf3, components=(30.,-20.))
bcs = (bc1, bc2, eLoad, nLoad, tLoad)

# nodes
# if one value is passed as parameter where oofem expects array of values, it must be passed as tuple or list (see load in n4)
n1 = oofemlib.node(1, domain, coords=(0., 0.,0. ), bc=(0,1,0))
n2 = oofemlib.node(2, domain, coords=(2.4,0.,0. ), bc=(0,0,0))
n3 = oofemlib.node(3, domain, coords=(3.8,0.,0. ), bc=(0,0,bc1))
n4 = oofemlib.node(4, domain, coords=(5.8,0.,1.5), bc=(0,0,0), load=(4,))
n5 = oofemlib.node(5, domain, coords=(7.8,0.,3.0), bc=(0,1,0))
n6 = oofemlib.node(6, domain, coords=(2.4,0.,3.0), bc=(bc1,1,bc2))
nodes = (n1, n2, n3, n4, n5, n6)

# material and cross section
mat = oofemlib.isoLE(1, domain, d=1., E=30.e6, n=0.2, tAlpha=1.2e-5)
cs  = oofemlib.simpleCS(1, domain, area=0.162, Iy=0.0039366, beamShearCoeff=1.e18, thick=0.54)

# elements
e1 = oofemlib.beam2d(1, domain, nodes=(1,n2),  mat=1,   crossSect=1,  boundaryLoads=(3,1), bodyLoads=(5,))
e2 = oofemlib.beam2d(2, domain, nodes=(2,3),   mat=mat, crossSect=1,  DofsToCondense=(6,), bodyLoads=[tLoad])
e3 = oofemlib.beam2d(3, domain, nodes=(n3,4),  mat=1,   crossSect=cs, dofstocondense=[3])
e4 = oofemlib.beam2d(4, domain, nodes=(n4,n5), mat=mat, crossSect=cs)
e5 = oofemlib.beam2d(5, domain, nodes=(n6,2),  mat=1,   crossSect=1,  DofsToCondense=(6,))
elems = (e1, e2, e3, e4, e5)

# add eveything to domain (resize container first to save some time, but it is not necessary 0 see ltfs)
domain.resizeDofManagers(len(nodes))
for n in nodes:
    domain.setDofManager(n.number, n)
domain.resizeElements(len(elems))
for e in elems:
    domain.setElement(e.number, e)
domain.setMaterial(1, mat)
domain.setCrossSection(1, cs)
domain.resizeBoundaryConditions(len(bcs))
for bc in bcs:
    domain.setBoundaryCondition(bc.number, bc)
for ltf in ltfs:
    domain.setLoadTimeFunction(ltf.number, ltf)


problem.checkProblemConsistency();
problem.setRenumberFlag();
problem.solveYourself();
problem.terminateAnalysis();
