#
# this example illustrates the user defined function in python injected into oofem
#
#
import oofempy

class MyFunc(oofempy.Function):
    def __init__(self, num, domain):
        oofempy.Function.__init__(self, num, domain)
    def evaluateAtTime (self, t):
        return 2.5*t
    def evaluateVelocityAtTime(self, t):
        return 0.0;
    def evaluateAccelerationAtTime(self, t):
        return 0.0;
    

# engngModel
problem = oofempy.linearStatic(nSteps=1, outFile="test3.out")

# domain (if no engngModel specified to domain, it is asigned to the last one created)
domain = oofempy.domain(1, 1, problem, oofempy.domainType._1dTrussMode, tstep_all=True, dofman_all=True, element_all=True)
problem.setDomain(1, domain, True)

ltf1 = oofempy.peakFunction(1, domain, t=1, f_t=1)
ltf2 = MyFunc(2, domain) # use custom ltf here
ltfs = (ltf1, ltf2)

# boundary conditions
# loadTimeFunction parameter can be specified as int value or as LoadTimeFunction itself (valid for all objects with giveNumber() method)
bc1   = oofempy.boundaryCondition(    1, domain, loadTimeFunction=1,    prescribedValue=0.0)
n2    = oofempy.nodalLoad(            2, domain, loadTimeFunction=2,    components=(1.,))
bcs = (bc1, n2)

# nodes
# if one value is passed as parameter where oofem expects array of values, it must be passed as tuple or list (see load in n4)
n1 = oofempy.node(1, domain, coords=(0.,  0., 0. ), bc=(bc1,))
n2 = oofempy.node(2, domain, coords=(2.4, 0., 0. ), load=(n2,))
nodes = (n1, n2)

# material and cross section
mat = oofempy.isoLE(1, domain, d=1., E=30.e6, n=0.2, tAlpha=1.2e-5)
cs  = oofempy.simpleCS(1, domain, area=0.5, Iy=0.0, beamShearCoeff=1.e18, thick=0.5)

# elements
e1 = oofempy.truss1d(1, domain, nodes=(1,n2),  mat=1,   crossSect=1)
elems = (e1,)

# add eveything to domain (resize container first to save some time, but it is not necessary 0 see ltfs)
domain.resizeDofManagers(len(nodes))
for n in nodes:
   domain.setDofManager(n.number, n)
domain.resizeElements(len(elems))
for e in elems:
   domain.setElement(e.number, e)
domain.resizeMaterials(1)
domain.setMaterial(1, mat)
domain.resizeCrossSectionModels(1)
domain.setCrossSection(1, cs)
domain.resizeBoundaryConditions(len(bcs))
for bc in bcs:
   domain.setBoundaryCondition(bc.number, bc)
domain.resizeFunctions(len(ltfs))
for ltf in ltfs:
   domain.setFunction(ltf.number, ltf)


print("\nSolving problem")
problem.checkProblemConsistency()
problem.init()
problem.postInitialize()
problem.setRenumberFlag()
problem.solveYourself()
problem.terminateAnalysis()
print("\nProblem solved")



