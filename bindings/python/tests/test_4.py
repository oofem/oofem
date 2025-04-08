#
# this example illustrates the user defined material model in python injected into oofem
#
#
try: # installed
    import oofem as oofempy
    from oofem import util
except: # in-tree
    import oofempy
    import util


class MyMaterial(oofempy.StructuralMaterial):
    def __init__(self, num, domain):
        oofempy.StructuralMaterial.__init__(self, num, domain)
        self.k = 1.5;
    def giveClassName(self):
        return "MyMaterial" 
    def giveInputRecordName(self):
        return "MyElement"
    # Overloading this method is not yet possible in pybind11
    # see https://github.com/pybind/pybind11/issues/1962
    # However a workaround is to override giveStatus, see below 
    def CreateStatus(self, gp):
        return oofempy.StructuralMaterialStatus (gp)
    def give1dStressStiffMtrx(self, mode, gp, tStep):
        answer = oofempy.FloatMatrix(1,1)
        answer[0,0] = self.k
        return answer;
    def giveRealStressVector_1d(self, reducedStrain, gp, tStep):
        answer = oofempy.FloatArray(1)
        answer[0] = self.k * reducedStrain[0]
        status = self.giveStatus(gp)
        status.letTempStrainVectorBe(reducedStrain);
        status.letTempStressVectorBe(answer);
        return answer;
    def giveStatus (self, gp):
        print ("getStatus")
        if (not gp.hasMaterialStatus()):
            print ("status created")
            status = oofempy.StructuralMaterialStatus (gp)
            print(status)
            gp.setMaterialStatus(status, oofempy.IntegrationPointStatusIDType.IPSID_Default)
            print(gp.giveMaterialStatus())
        return gp.giveMaterialStatus() 
        
def test_4():
    # engngModel
    problem = oofempy.linearStatic(nSteps=1, outFile="test_4.out")

    # domain (if no engngModel specified to domain, it is asigned to the last one created)
    domain = oofempy.domain(1, 1, problem, oofempy.domainType._1dTrussMode, tstep_all=True, dofman_all=True, element_all=True)
    problem.setDomain(1, domain, True)

    ltf1 = oofempy.peakFunction(1, domain, t=1, f_t=1)
    ltfs = (ltf1, )

    # boundary conditions
    # loadTimeFunction parameter can be specified as int value or as LoadTimeFunction itself (valid for all objects with giveNumber() method)
    bc1   = oofempy.boundaryCondition(    1, domain, loadTimeFunction=1,    prescribedValue=0.0)
    n2    = oofempy.nodalLoad(            2, domain, loadTimeFunction=1,    components=(1.,), dofs=(1,))
    bcs = (bc1, n2)

    # nodes
    # if one value is passed as parameter where oofem expects array of values, it must be passed as tuple or list (see load in n4)
    n1 = oofempy.node(1, domain, coords=(0.,  0., 0. ), bc=(bc1,))
    n2 = oofempy.node(2, domain, coords=(2.4, 0., 0. ), load=(n2,))
    nodes = (n1, n2)

    # material and cross section
    #mat = oofempy.isoLE(1, domain, d=1., E=30.e6, n=0.2, tAlpha=1.2e-5)
    mat = MyMaterial(1, domain)
    cs  = oofempy.simpleCS(1, domain, area=0.5, Iy=0.0, beamShearCoeff=1.e18, thick=0.5)

    # elements
    e1 = oofempy.truss1d(1, domain, nodes=(1,n2),  mat=1,   crossSect=1)
    # construct OOFEMTXTInputRecord from bp::dict **kw
    elems = (e1,)

    # setup domain
    util.setupDomain (domain, nodes, elems, (cs,), (mat,), bcs, (), ltfs, ())

    print("\nSolving problem")
    problem.checkProblemConsistency()
    problem.init()
    problem.postInitialize()

    # Test user defined material model 
    a =oofempy.FloatArray((1.0,))
    ans=oofempy.FloatArray()

    rule=e1.giveDefaultIntegrationRulePtr()
    rule=e1.giveIntegrationRule(0)
    #print(rule)
    gp = rule.getIntegrationPoint(0)
    #print(gp)
    print("Printing status", mat.giveStatus(gp))
    ans = domain.giveMaterial(1).giveRealStressVector_1d(a, rule.getIntegrationPoint(0), None)
    #ans.pY()
    assert (round(ans[0]-1.5, 8) == 0), "Stress value error"

    problem.setRenumberFlag()
    problem.solveYourself()

    #check solution
    u2 = problem.giveUnknownComponent (oofempy.ValueModeType.VM_Total, problem.giveCurrentStep(False), domain, domain.giveDofManager(2).giveDofWithID(oofempy.DofIDItem.D_u))
    assert (round (u2-3.2, 8) == 0), "Node 2 dof 1 displacement check failed"


    problem.terminateAnalysis()
    print("\nProblem solved")



if __name__ == "__main__":
    test_4()

