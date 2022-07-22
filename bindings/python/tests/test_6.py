# this example illustrates an injection of external temperature field during runtime
import oofempy
import util


def test_6():
    # engngModel
    problem = oofempy.staticStructural(nSteps=1, outFile="test_6.out")

    # domain (if no engngModel specified to domain, it is asigned to the last one created)
    domain = oofempy.domain(1, 1, problem, oofempy.domainType._3dMode, tstep_all=True, dofman_all=True, element_all=True)
    problem.setDomain(1, domain, True)
    # load time functions
    ltf1 = oofempy.constantFunction(1, domain, f_t=1.)
    ltfs = (ltf1,)
    # boundary conditions
    bc1 = oofempy.boundaryCondition(1, domain, loadTimeFunction=1, prescribedValue=0.0)
    dw = oofempy.DeadWeight(2, domain, loadTimeFunction = 1, components=(0.,0.,-10.0))
    bcs = (bc1, dw)
    # nodes
    n1 = oofempy.node(1, domain, coords=(0., 0., 2. ))
    n2 = oofempy.node(2, domain, coords=(0., 2., 2. ))
    n3 = oofempy.node(3, domain, coords=(2., 2., 2. ))
    n4 = oofempy.node(4, domain, coords=(2., 0., 2. ))
    n5 = oofempy.node(5, domain, coords=(0., 0., 0. ), bc=(bc1, bc1, bc1))
    n6 = oofempy.node(6, domain, coords=(0., 2., 0. ), bc=(  0,   0, bc1) )
    n7 = oofempy.node(7, domain, coords=(2., 2., 0. ), bc=(  0,   0, bc1))
    n8 = oofempy.node(8, domain, coords=(2., 0., 0. ), bc=(  0, bc1, bc1))
    #
    nodes = (n1,n2,n3,n4,n5,n6,n7,n8)
    # material and cross section
    mat = oofempy.isoLE(1, domain, d=2., E=100., n=0.01, tAlpha=1.2e-5)
    cs  = oofempy.simpleCS(1, domain)
    # elements - assign through sets
    e1 = oofempy.lspace(1, domain, nodes=(1,2,3,4,5,6,7,8), mat=1, crossSect = 1, bodyloads = (dw,))
    elems = (e1,)
    #sets
    set1 = oofempy.createSet(1, domain, elements=(1,))
    # setup domain
    util.setupDomain(domain, nodes, elems, (cs,), (mat,), bcs, (), ltfs, (set1,))
    # add export module for outputting regular VTU files - automatically registered
#    vtkxml = oofempy.vtkxml(1, problem, domain_all=True, tstep_all=True, dofman_all=True, element_all=True, vars=(56,), primvars=(6,), stype=1, pythonExport=0)
   
    # add export module for outputting python lists with values - automatically registered
    print("\nSolving problem")
    problem.checkProblemConsistency()
    problem.init()
    problem.postInitialize()
    problem.giveTimer().startTimer(oofempy.EngngModelTimerType.EMTT_AnalysisTimer)
    activeMStep = problem.giveMetaStep(1)
    problem.initMetaStepAttributes(activeMStep);
       
    for timeStep in range(1):
        problem.preInitializeNextStep()
        problem.giveNextStep()
        currentStep = problem.giveCurrentStep()
        problem.initializeYourself( currentStep )
        problem.solveYourselfAt( currentStep )
        problem.updateYourself( currentStep )
        problem.terminate( currentStep )
        print("TimeStep %d finished" % (timeStep))

    #check solution
    w3 = problem.giveUnknownComponent(oofempy.ValueModeType.VM_Total, problem.giveCurrentStep(False), domain, domain.giveDofManager(1).giveDofWithID(oofempy.DofIDItem.D_w))
    assert (round (w3+4.000e-1, 8) == 0), "Node 1 dof 3 displacement check failed"
    
        
    problem.terminateAnalysis()
    print("\nProblem solved")


if __name__ == "__main__":
    test_6()
