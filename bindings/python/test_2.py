import oofempy
import util

def test_2():

    # engngModel
    problem = oofempy.linearStatic(nSteps=3, outFile='test_2.out')

    # domain (if no engngModel specified to domain, it is asigned to the last one created)
    domain = oofempy.domain(1, 1, problem, oofempy.domainType._2dBeamMode, tstep_all=True, dofman_all=True, element_all=True)
    problem.setDomain(1, domain, True)

    ltf1 = oofempy.peakFunction(1, domain, t=1, f_t=1)
    ltf2 = oofempy.peakFunction(2, domain, t=2, f_t=1)
    ltf3 = oofempy.peakFunction(3, domain, t=3, f_t=1)
    ltfs = (ltf1, ltf2, ltf3)

    # boundary conditions
    # loadTimeFunction parameter can be specified as int value or as LoadTimeFunction itself (valid for all objects with giveNumber() method)
    bc1   = oofempy.boundaryCondition(    1, domain, loadTimeFunction=1,    prescribedValue=0.0)
    bc2   = oofempy.boundaryCondition(    2, domain, loadTimeFunction=2,    prescribedValue=-.006e-3)
    eLoad = oofempy.constantEdgeLoad(     3, domain, loadTimeFunction=1, components=(0.,10.,0.), loadType=3, ndofs=3)
    nLoad = oofempy.nodalLoad(            4, domain, loadTimeFunction=1,    components=(-18.,24.,0.))
    tLoad = oofempy.structTemperatureLoad(5, domain, loadTimeFunction=3, components=(30.,-20.))
    bcs = (bc1, bc2, eLoad, nLoad, tLoad)

    # nodes
    # if one value is passed as parameter where oofem expects array of values, it must be passed as tuple or list (see load in n4)
    n1 = oofempy.node(1, domain, coords=(0.,  0., 0. ), bc=(0,1,0))
    n2 = oofempy.node(2, domain, coords=(2.4, 0., 0. ), bc=(0,0,0))
    n3 = oofempy.node(3, domain, coords=(3.8, 0., 0. ), bc=(0,0,bc1))
    n4 = oofempy.node(4, domain, coords=(5.8, 0., 1.5), bc=(0,0,0), load=(4,))
    n5 = oofempy.node(5, domain, coords=(7.8, 0., 3.0), bc=(0,1,0))
    n6 = oofempy.node(6, domain, coords=(2.4, 0., 3.0), bc=(bc1,1,bc2))
    nodes = (n1, n2, n3, n4, n5, n6)

    # material and cross section
    mat = oofempy.isoLE(1, domain, d=1., E=30.e6, n=0.2, tAlpha=1.2e-5)
    cs  = oofempy.simpleCS(1, domain, area=0.162, Iy=0.0039366, beamShearCoeff=1.e18, thick=0.54)

    # elements
    e1 = oofempy.beam2d(1, domain, nodes=(1,n2),  mat=1,   crossSect=1,  boundaryLoads=(3,1), bodyLoads=(5,))
    e2 = oofempy.beam2d(2, domain, nodes=(2,3),   mat=mat, crossSect=1,  DofsToCondense=(6,), bodyLoads=[tLoad])
    e3 = oofempy.beam2d(3, domain, nodes=(n3,4),  mat=1,   crossSect=cs, dofstocondense=[3])
    e4 = oofempy.beam2d(4, domain, nodes=(n4,n5), mat=mat, crossSect=cs)
    e5 = oofempy.beam2d(5, domain, nodes=(n6,2),  mat=1,   crossSect=1,  DofsToCondense=(6,))
    elems = (e1, e2, e3, e4, e5)

    # add eveything to domain (resize container first to save some time, but it is not necessary 0 see ltfs)
    util.setupDomain(domain, nodes, elems, (cs,), (mat,),  bcs, (), ltfs, ())

    print("\nSolving problem")
    problem.checkProblemConsistency()
    problem.init()
    problem.postInitialize()
    problem.setRenumberFlag()
    problem.solveYourself()

    #check solution
    u1 = problem.giveUnknownComponent (oofempy.ValueModeType.VM_Total, problem.giveCurrentStep(False), domain, domain.giveDofManager(1).giveDofWithID(oofempy.DofIDItem.D_u))
    print (u1)
    assert (round (u1+8.64000000e-04, 8) == 0), "Node 1 dof 1 displacement check failed"
    u4 = problem.giveUnknownComponent (oofempy.ValueModeType.VM_Total, problem.giveCurrentStep(False), domain, domain.giveDofManager(4).giveDofWithID(oofempy.DofIDItem.D_u))
    assert (round (u4-9.47333333e-04, 8) == 0), "Node 4 dof 1 displacement check failed"
    

    problem.terminateAnalysis()
    print("\nProblem solved")

if __name__ == "__main__":
    test_2()
