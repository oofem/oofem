#
# this example illustrates the user defined function in python injected into oofem
#
#
try: # installed
    import oofem as oofempy
    from oofem import util
except: # in-tree
    import oofempy
    import util
import numpy as np


class MyFunc(oofempy.Function):
    def __init__(self, num, domain):
        oofempy.Function.__init__(self, num, domain)
    def evaluateAtTime (self, t):
        return 2.5*t
    def evaluateVelocityAtTime(self, t):
        return 0.0;
    def evaluateAccelerationAtTime(self, t):
        return 0.0;

class MyElement(oofempy.StructuralElement):
    def __init__(self, num, domain):
        oofempy.StructuralElement.__init__(self, num, domain)
        self.setNumberOfDofManagers(2)
    def computeStiffnessMatrix(self, answer, rMode, tStep):
        answer.resize(2,2);
        answer[0,0] = 1.0;
        answer[0,1] = -1.0;
        answer[1,0] = 1.0;
        answer[1,1] = -1.0;
    def computeBmatrixAt(self, gp, answer, lowerIndx, upperIndx):
        answer.resize(1,2);
        answer[0,0]=-1.0;
        answer[0,1]=1.0;
    def giveInternalForcesVector(self, answer, tStep, useUpdatedGpRecord):
        u = oofempy.FloatArray()
        k = oofempy.FloatMatrix(2,2)
        self.computeVectorOf(oofempy.ValueModeType.VM_Total, tStep, u)
        self.computeStiffnessMatrix(k, oofempy.CharType.StiffnessMatrix, tStep)
        answer.beProductOf (k, u)
    def giveNumberOfDofs(self):
        return 2;
    def computeNumberOfDofs(self):
        return 2;
    def giveDofManDofIDMask(self, inode, answer):
        answer.resize(1)
        answer[0] = oofempy.DofIDItem.D_u
#    def printOutputAt(self, f, tStep):
#        fint = oofempy.FloatArray()
#        self.giveInternalForcesVector(fint, tStep, false)
#        fint.py()
    def giveClassName(self):
        return "MyElement"
    def giveInputRecordName(self):
        return "MyElement"

def test_3():
    # engngModel
    problem = oofempy.linearStatic(nSteps=1, outFile="test_3.out")

    # domain (if no engngModel specified to domain, it is asigned to the last one created)
    domain = oofempy.domain(1, 1, problem, oofempy.domainType._1dTrussMode, tstep_all=True, dofman_all=True, element_all=True)
    problem.setDomain(1, domain, True)

    ltf1 = oofempy.peakFunction(1, domain, t=1, f_t=1)
    ltf2 = MyFunc(2, domain) # use custom ltf here
    ltfs = (ltf1, ltf2)

    # boundary conditions
    # loadTimeFunction parameter can be specified as int value or as LoadTimeFunction itself (valid for all objects with giveNumber() method)
    bc1   = oofempy.boundaryCondition(    1, domain, loadTimeFunction=1,    prescribedValue=0.0)
    n2    = oofempy.nodalLoad(            2, domain, loadTimeFunction=2,    components=(1.,), dofs=(1,))
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
    #e1 = oofempy.truss1d(1, domain, nodes=(1,n2),  mat=1,   crossSect=1)
    e1 = MyElement(1, domain)
    e1.setDofManagers(np.array([1,2]))
    # construct OOFEMTXTInputRecord from bp::dict **kw
    ir = oofempy.OOFEMTXTInputRecord()
    ir.setRecordString ("nodes 2 1 2 mat 1 crosssect 1")
    # pass input record to elem
    e1.initializeFrom(ir,2)
    elems = (e1,)

    # setup domain
    util.setupDomain (domain, nodes, elems, (cs,), (mat,), bcs, (), ltfs, ())

    print("\nInitializing problem")
    problem.checkProblemConsistency()
    problem.init()
    problem.postInitialize()
    print("\nSolving problem")
    #problem.setRenumberFlag()
    problem.solveYourself()
    print("\nProblem solved")

    #check solution
    u2 = problem.giveUnknownComponent (oofempy.ValueModeType.VM_Total, problem.giveCurrentStep(False), domain, domain.giveDofManager(2).giveDofWithID(oofempy.DofIDItem.D_u))
    assert (round (u2+2.5, 8) == 0), "Node 2 dof 1 displacement check failed"

    problem.terminateAnalysis()
    print("\nProblem solved")


if __name__ == "__main__":
    test_3()
