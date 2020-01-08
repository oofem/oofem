#
# this example illustrates the user defined material model in python injected into oofem
#
#
from __future__ import print_function
import numpy as np
import sys
sys.path.insert(0, '../../default')
import oofempy
import util


class MyMaterial(oofempy.StructuralMaterial):
    def __init__(self, num, domain):
        oofempy.StructuralMaterial.__init__(self, num, domain)
        self.k = 1.5
        self.talpha = 0.00001
        self.referencetemperature = 0.

    def giveClassName(self):
        return "MyMaterial"

    def giveInputRecordName(self):
        return "MyElement"

    # Overloading this method is not yet possible in pybind11
    # see https://github.com/pybind/pybind11/issues/1962
    # However a workaround is to override giveStatus, see below 
    def CreateStatus(self, gp):
        return oofempy.StructuralMaterialStatus(gp)

    def give1dStressStiffMtrx(self, mode, gp, tStep):
        answer = oofempy.FloatMatrix(1, 1)
        answer[0, 0] = self.k
        return answer

    def giveRealStressVector_1d(self, reducedStrain, gp, tStep):
        answer = oofempy.FloatArray(1)
        answer[0] = self.k * reducedStrain[0]
        status = self.giveStatus(gp)
        status.letTempStrainVectorBe(reducedStrain)
        status.letTempStressVectorBe(answer)
        return answer

    def giveStatus(self, gp):
        print("getStatus")
        if gp.giveMaterialStatus() is None:
            print("getStatus creating")
            status = oofempy.StructuralMaterialStatus(gp)
            gp.setMaterialStatus(status)
        return gp.giveMaterialStatus() 


def test_5():

    f = oofempy.UniformGridField()
    dr = oofempy.OOFEMTXTDataReader('input_5.in')
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
