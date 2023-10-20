import oofempy
import util

#Transient transport problem with one asisymmetric element loaded by heat flux from the right edge

def test_7():
    # engngModel
    problem = oofempy.transientTransport(nSteps=10, deltat=3600.0, alpha=0.5, outFile='test_7.out')
    domain = oofempy.domain(1, 1, problem, oofempy.domainType._HeatTransferMode, tstep_all=True, dofman_all=True, element_all=True)
    problem.setDomain(1, domain, True)
   
    #nodes
    n1 = oofempy.node(1, domain, coords=(0.1, 0.3, 0. ))
    n2 = oofempy.node(2, domain, coords=(1.1, 0.3, 0. ))
    n3 = oofempy.node(3, domain, coords=(1.1, 1.5, 0. ))
    n4 = oofempy.node(4, domain, coords=(0.1, 1.5, 0. ))
   
    # elements
    e1 = oofempy.quadAxiSym1ht(1, domain, nodes=(1,2,3,4), mat=1, crossSect=1)
   
    # material and cross section
    mat = oofempy.isoHeat(1, domain, d=2400., k=1.5, c=800.)
    cs  = oofempy.simpleTransportCS(1, domain, mat=1, thickness=0.1)

    # boundary conditions
    bc1 = oofempy.boundaryCondition(1, domain, loadTimeFunction=1, dofs=(10,), values=(50.,), set=1)
    bc2 = oofempy.constantEdgeLoad(2, domain, loadTimeFunction=1, components=(-150.0,), loadtype=2,set=2)

    #Initial condition
    ic1 = oofempy.initialCondition(1, domain, Conditions=1, u=50.0, dofs=(10,), set=3)

    #Load time functions
    ltf1 = oofempy.constantFunction(1, domain, f_t=1.)

    #Sets
    set1 = oofempy.createSet(1, domain, nodes=(1,4))
    set2 = oofempy.createSet(2, domain, elementedges=(1,2))
    set3 = oofempy.createSet(3, domain, nodes=(1,2,3,4))

    # add eveything to domain
    util.setupDomain(domain, nodes=(n1,n2,n3,n4), elems=(e1,), css=(cs,), mats=(mat,), bcs=(bc1,bc2), ics=(ic1,), ltfs=(ltf1,), sets=(set1,set2,set3))

    problem.checkProblemConsistency()
    problem.init()
    problem.postInitialize()

    problem.giveTimer().startTimer(oofempy.EngngModelTimerType.EMTT_AnalysisTimer)
    activeMStep = problem.giveMetaStep(1)
    problem.initMetaStepAttributes(activeMStep);
    problem.setRenumberFlag()
    #vtkxml = oofempy.vtkxml(1, problem, domain_all=True, tstep_all=True, dofman_all=True, element_all=True, primvars=(6,), stype=1, pythonExport=0)
    for step in range(5): 
        problem.preInitializeNextStep()
        problem.giveNextStep()
        problem.forceEquationNumbering()
        currentStep = problem.giveCurrentStep()
        currentStep.setTargetTime(float(step))
        problem.initializeYourself( currentStep )
        problem.solveYourselfAt( currentStep )
        problem.updateYourself( currentStep )
        problem.terminate( currentStep )
    #check solution
    t1 = problem.giveUnknownComponent(oofempy.ValueModeType.VM_Total, problem.giveCurrentStep(False), domain, domain.giveDofManager(2).giveDofWithID(10))
    assert (round (t1-55.37908, 4) == 0), "Temperature at node 2, check failed"
    
    problem.terminateAnalysis()
    print("\nProblem solved")

if __name__ == "__main__":
    test_7()


