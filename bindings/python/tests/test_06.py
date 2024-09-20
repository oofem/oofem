# this example illustrates an injection of external eigenstrain
import oofempy
import util

def test_06():
    # engngModel
    problem = oofempy.staticStructural(nSteps=1, outFile="test_06.out")

    # domain (if no engngModel specified to domain, it is asigned to the last one created)
    domain = oofempy.domain(1, 1, problem, oofempy.domainType._2dPlaneStressMode, tstep_all=True, dofman_all=True, element_all=True)
    problem.setDomain(1, domain, True)
    
    ltf1 = oofempy.constantFunction(1, domain, f_t=1.)
    ltf2 = oofempy.piecewiseLinFunction(2, domain, t=(1., 5.), f_t= (0., 4.))
    ltfs = (ltf1,ltf2)

    # boundary conditions
    bc1 = oofempy.boundaryCondition(1, domain, loadTimeFunction=1, prescribedValue=0.0)
    bcs = (bc1,)

    # nodes
    n1 = oofempy.node(1, domain, coords=(0.,  0.4, 0. ),)
    n2 = oofempy.node(2, domain, coords=(2.4, 0., 0. ), bc=(bc1,bc1))
    n3 = oofempy.node(3, domain, coords=(2.4, 0.8, 0.), bc=(bc1,bc1))
    n4 = oofempy.node(4, domain, coords=(4.8, 0.4, 0.), )
    nodes = (n1,n2,n3,n4)
    
    # material and cross section
    mat = oofempy.isoLE(1, domain, d=1., E=30.e3, n=0.2, tAlpha=0.)
    cs  = oofempy.simpleCS(1, domain, thick=0.5, material=1, set=1)
    
    # elements - assign through sets
    e1 = oofempy.trPlaneStress2d(1, domain, nodes=(1,2,3), nip=4)
    e2 = oofempy.trPlaneStress2d(2, domain, nodes=(2,4,3), nip=4)
    elems = (e1,e2)

    set1 = oofempy.createSet(1, domain, elements=(1,2))
    set2 = oofempy.createSet(2, domain, elements=(1,))
    sets = (set1, set2)

    # setup domain
    util.setupDomain(domain, nodes, elems, (cs,), (mat,), bcs, (), ltfs, sets)
    
    # add export module for outputting regular VTU files - automatically registered
    vtkxml = oofempy.vtkxml(1, problem, domain_all=True, tstep_all=True, dofman_all=True, element_all=True, vars=(1,4), primvars=(1,), stype=1, pythonExport=0)
    
    # add export module for outputting python lists with values - automatically registered
    vtkPy = oofempy.vtkmemory(1, problem, domain_all=True, tstep_all=True, dofman_all=True, element_all=True, vars=(1,4), primvars=(1,), cellvars = (47,103), stype=1, pythonExport=1)
    
        
    print("\nSolving problem")
    problem.checkProblemConsistency()
    problem.init()
    problem.postInitialize()
    #print(e1.giveMaterial())
    problem.giveTimer().startTimer(oofempy.EngngModelTimerType.EMTT_AnalysisTimer)
    activeMStep = problem.giveMetaStep(1)
    problem.initMetaStepAttributes(activeMStep);
    
    #Create dummy eigenstrain field
    f = oofempy.UniformGridField()
    f.setGeometry((0., 0., 0.), (6., 1., 1), (1, 1, 1)) #(low), (high), (division)
    vals = []
    for x in range (0,2):
        for y in range (0,2):
            for z in range (0,2):
                vals.append(oofempy.FloatArray((1.e-3*x, 2.e-4, 0.)))
    f.setValues(vals) #Eigenstrain, linear in x direction
    f.setSetsNumbers((2,)) #Assign field only to element "1" defined in set "2"
    context = problem.giveContext()
    field_man = context.giveFieldManager()
    field_man.registerField(f, oofempy.FieldType.FT_EigenStrain)
    
    for timeStep in range(1):
        problem.preInitializeNextStep()
        problem.giveNextStep()
        currentStep = problem.giveCurrentStep()
        problem.initializeYourself( currentStep )
        problem.solveYourselfAt( currentStep )
        problem.updateYourself( currentStep )
        problem.terminate( currentStep )
        print("TimeStep %d finished" % (timeStep))
        answer = oofempy.FloatArray()
        coord = (3.0, 0., 0.9)
        f.evaluateAt(answer, coord, oofempy.ValueModeType.VM_Total, currentStep)
        print("Interpolated eigenstrain at", coord, answer)

        #Nodal values, obtained by smoothing to nodes (defined by stype=1). Discontinuous fields (e.g. stress) are made continous.
        for p in vtkPy.getExportRegions():
            print ("Piece:", p)
            print("Vertices:", p.getVertices())
            print("Cells:", p.getCellConnectivity())
            print("CellTypes:", p.getCellTypes())
            displ = p.getPrimaryVertexValues(oofempy.UnknownType.DisplacementVector);
            print ("Displacement in vertices:", displ)
            d1=p.givePrimaryVarInNode(oofempy.UnknownType.DisplacementVector,1)
            print ("Displacement in node 1:", d1)
            assert (round (d1[0]+0.000736, 8) == 0), "Node 1 dof 1 displacement check failed, should be -0.000736"
            d3=p.givePrimaryVarInNode(oofempy.UnknownType.DisplacementVector,3)
            print ("Displacement in node 3:", d3)
            assert (round (d3[0], 8) == 0), "Node 3 dof 1 displacement check failed, should be 0.0"
            print ("Stress in node 1:", p.giveInternalVarInNode(oofempy.InternalStateType.IST_StressTensor,1))
            print ("Strain in node 1:", p.giveInternalVarInNode(oofempy.InternalStateType.IST_StrainTensor,1))

        #Extract field from a nearest integration point. No interpolation.
        sl = domain.giveSpatialLocalizer()
        coord = oofempy.FloatArray((1.2, 0.5, 0))
        ip = sl.giveClosestIP(coord, set1, 0)
        IPgCoords = ip.giveGlobalCoordinates()
        elem = ip.giveElement()
        mat=elem.giveCrossSection().giveMaterial(ip)
        elem.giveIPValue(answer,ip, oofempy.InternalStateType.IST_StressTensor, currentStep)
        print("Stress in IP at:", IPgCoords, answer)
        elem.giveIPValue(answer,ip, oofempy.InternalStateType.IST_StrainTensor, currentStep)
        print("Strain in IP at:", IPgCoords, answer)

        #Set internal variable field and interpolate on element using shape functions.
        internalFieldStress = oofempy.InternalVariableField( oofempy.InternalStateType.IST_StressTensor, oofempy.FieldType.FT_Unknown, oofempy.MaterialMappingAlgorithmType.MMA_ShapeFunctionProjection, domain)   
        answer=oofempy.FloatArray()
        internalFieldStress.evaluateAt(answer, coord, oofempy.ValueModeType.VM_Total, currentStep)
        print("Stress from InternalVariableField at:", coord, answer)


    problem.terminateAnalysis()
    

if __name__ == "__main__":
    test_06()
