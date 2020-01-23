import oofempy
#Inject external temperature field to OOFEM problem, defined on a uniform grid

def example_1():
    dr = oofempy.OOFEMTXTDataReader('example_1.in')
    problem = oofempy.InstanciateProblem(dr, oofempy.problemMode.processor, 0, None, False)

    f = oofempy.UniformGridField()    
    #Lower bound, upper bound, number of elements on edges
    f.setGeometry((0., 0., 0.), (1., 1., 1.), (2, 2, 2))
    f.setValues([500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500.])

    context = problem.giveContext()
    field_man = context.giveFieldManager()
    field_man.registerField(f, oofempy.FieldType.FT_Temperature)

    print("\nSolving problem")
    problem.checkProblemConsistency()
    problem.init()
    problem.postInitialize()

    problem.setRenumberFlag()
    problem.solveYourself()

    problem.terminateAnalysis()
    print("\nProblem solved")


if __name__ == "__main__":
    example_1()
