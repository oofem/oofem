import oofempy
#Inject external temperature field to OOFEM problem, defined on a uniform grid

def example_1():
    dr = oofempy.OOFEMTXTDataReader('example_1.in')
    problem = oofempy.InstanciateProblem(dr, oofempy.problemMode.processor, 0, None, False)

    # add export module for outputting python lists with values - automatically registered
    vtkxmlPy = oofempy.vtkxml(1, problem, domain_all=True, tstep_all=True, dofman_all=True, element_all=True, vars=(1,4), primvars=(1,), stype=1, pythonExport=1)

    f = oofempy.UniformGridField()    
    #Lower bound, upper bound, number of elements on edges
    f.setGeometry((0., 0., 0.), (1., 1., 1.), (2, 2, 2))
    f.setValues(list(range(27))) #put some temperature in the injected field

    context = problem.giveContext()
    field_man = context.giveFieldManager()
    field_man.registerField(f, oofempy.FieldType.FT_Temperature)

    print("\nSolving problem")
    problem.checkProblemConsistency()
    problem.init()
    problem.postInitialize()

    problem.setRenumberFlag()
    problem.solveYourself()

    print("Nodes", vtkxmlPy.getNodes()['1'])
    print("Displacements", vtkxmlPy.getPrimaryVars()['DisplacementVector'])    
    print("Stress tensor at node 1:", vtkxmlPy.getInternalVars()['IST_StressTensor'][0])
    print("Strain tensor at node 1:", vtkxmlPy.getInternalVars()['IST_StrainTensor'][0])


    problem.terminateAnalysis()
    print("\nProblem solved")
   

if __name__ == "__main__":
    example_1()
