from __future__ import print_function
import oofempy


def test_5():

    f = oofempy.UniformGridField()
    dr = oofempy.OOFEMTXTDataReader('example_1.in')
    problem = oofempy.InstanciateProblem(dr, oofempy.problemMode.processor, 0, None, False)

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
    test_5()
