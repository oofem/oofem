/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <Python.h>

#include "structuralpythonmaterial.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(StructuralPythonMaterial);

StructuralPythonMaterial :: StructuralPythonMaterial(int n, Domain *d) :
    StructuralMaterial(n, d)
{}

StructuralPythonMaterial :: ~StructuralPythonMaterial()
{
///    if ( mpModule ) Py_DECREF(mpModule);
}

IRResultType StructuralPythonMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, this->moduleName, _IFT_StructuralPythonMaterial_moduleName);

    module=bp::import(moduleName.c_str());
    if(!module){ OOFEM_WARNING("Module %s not importable.",moduleName.c_str()); return IRRT_BAD_FORMAT; }
    // lambda for finding function and checking that it is callable
    // returns true (OK) if the function was not found, or was found and it callable
    // returns false (not OK) if the function was found but is not callable
    auto tryDef=[&](const std::string& attr, bp::object& callable)->bool{
        if(PyObject_HasAttrString(module.ptr(),attr.c_str())){
            callable=module.attr(attr.c_str());
            if(!PyCallable_Check(callable.ptr())){ OOFEM_WARNING("Object %s is not callable",attr.c_str()); return false; }
        } else callable=bp::object();
        return true;
    };
    // try to find all necessary functions; false means the function is not callable, in which case warning was already printed above
    if(!(tryDef("computeStress",smallDef) && tryDef("computePK1Stress",largeDef) && tryDef("computeStressTangent",smallDefTangent) && tryDef("computePK1StressTangent",largeDefTangent))){ return IRRT_BAD_FORMAT; }
    if(!smallDefTangent && !!smallDef){ OOFEM_WARNING("Using numerical tangent for small deformations."); }
    if(!largeDefTangent && !!largeDef){ OOFEM_WARNING("Using numerical tangent for large deformations."); }
    if(!smallDef && !largeDef){ OOFEM_WARNING("No functions for small/large deformations found."); return IRRT_BAD_FORMAT; }

    pert = 1e-12;

    return IRRT_OK;
}

void StructuralPythonMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);

    input.setField(this->moduleName, _IFT_StructuralPythonMaterial_moduleName);
}

MaterialStatus *StructuralPythonMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralPythonMaterialStatus(gp);
}

void StructuralPythonMaterial :: callStressFunction(bp::object func, const FloatArray &oldStrain, const FloatArray &oldStress, const FloatArray &strain, FloatArray &stress, bp::object stateDict, bp::object tempStateDict, TimeStep *tStep) const
{
    // pass mutable args via bp::ref
    // pass "const" args without, which by default results in a new copy, ensuring the original won't be modified
    func(oldStrain,oldStress,strain,stress,stateDict,tempStateDict,tStep->giveTargetTime());
}

void StructuralPythonMaterial :: callTangentFunction(FloatMatrix &answer, bp::object func, const FloatArray &strain, const FloatArray &stress, bp::object stateDict, bp::object tempStateDict, TimeStep *tStep) const
{
    answer=bp::extract<FloatMatrix>(func(strain,stress,stateDict,tempStateDict,tStep->giveTargetTime()));
}

void StructuralPythonMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralPythonMaterialStatus *ms = dynamic_cast< StructuralPythonMaterialStatus * >( this->giveStatus(gp) );

    if ( this->smallDefTangent ) {
        this->callTangentFunction(answer, this->smallDefTangent, ms->giveTempStrainVector(), ms->giveTempStressVector(), ms->giveStateDictionary(), ms->giveTempStateDictionary(), tStep);
    } else {
        FloatArray vE, vE_h, stress, stressh;
        vE = ms->giveTempStrainVector();
        stress = ( ( StructuralMaterialStatus * ) gp->giveMaterialStatus() )->giveTempStressVector();
        answer.resize(6, 6);
        for ( int i = 1; i <= 6; ++i ) {
            vE_h = vE;
            vE_h.at(i) += pert;
            this->giveRealStressVector_3d(stressh, gp, vE_h, tStep);
            stressh.subtract(stress);
            stressh.times(1.0 / pert);
            answer.setColumn(stressh, i);
        }

        // Reset the stress internal variables
        this->giveRealStressVector_3d(stress, gp, vE, tStep);
    }
}


void StructuralPythonMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralPythonMaterialStatus *ms = dynamic_cast< StructuralPythonMaterialStatus * >( this->giveStatus(gp) );

    if ( this->largeDefTangent ) {
        this->callTangentFunction(answer, this->largeDefTangent, ms->giveTempFVector(), ms->giveTempPVector(), ms->giveStateDictionary(), ms->giveTempStateDictionary(), tStep);
    } else {
        FloatArray vF, vF_h, stress, stressh;
        vF = ms->giveTempFVector();
        stress = ( ( StructuralMaterialStatus * ) gp->giveMaterialStatus() )->giveTempPVector();
        answer.resize(9, 9);
        for ( int i = 1; i <= 9; ++i ) {
            vF_h = vF;
            vF_h.at(i) += pert;
            this->giveFirstPKStressVector_3d(stressh, gp, vF_h, tStep);
            stressh.subtract(stress);
            stressh.times(1.0 / pert);
            answer.setColumn(stressh, i);
        }

        // Reset the internal variables
        this->giveFirstPKStressVector_3d(stress, gp, vF, tStep);
    }
}


void StructuralPythonMaterial :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                                   const FloatArray &strain, TimeStep *tStep)
{
    StructuralPythonMaterialStatus *ms = static_cast< StructuralPythonMaterialStatus * >( this->giveStatus(gp) );

    ms->reinitTempStateDictionary();

    this->callStressFunction(this->smallDef, 
                              ms->giveStrainVector(), ms->giveStressVector(),
                              strain, answer,
                              ms->giveStateDictionary(), ms->giveTempStateDictionary(), tStep);

    ms->letTempStrainVectorBe(strain);
    ms->letTempStressVectorBe(answer);
}


void StructuralPythonMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                                      const FloatArray &vF, TimeStep *tStep)
{
    StructuralPythonMaterialStatus *ms = static_cast< StructuralPythonMaterialStatus * >( this->giveStatus(gp) );

    ms->reinitTempStateDictionary(); // Resets the temp dictionary to the equilibrated values

    this->callStressFunction(this->smallDef, 
                            ms->giveFVector(), ms->givePVector(),
                            vF, answer,
                            ms->giveStateDictionary(), ms->giveTempStateDictionary(), tStep);

    FloatArray vE, vS;
    FloatMatrix F, Finv, E, S;
    F.beMatrixForm(vF);
    Finv.beInverseOf(F);
    // Compute Green-Lagrange strain
    E.beTProductOf(F, F);
    E.at(1, 1) -= 1.0;
    E.at(2, 2) -= 1.0;
    E.at(3, 3) -= 1.0;
    E.times(0.5);
    vE.beSymVectorFormOfStrain(E);
    // Convert from P to S
    S.beProductOf(Finv, answer);
    vS.beSymVectorForm(S);

    ms->letTempStrainVectorBe(vE);
    ms->letTempStressVectorBe(vS);
    ms->letTempPVectorBe(answer);
    ms->letTempFVectorBe(vF);
}


int StructuralPythonMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    StructuralPythonMaterialStatus *ms = static_cast< StructuralPythonMaterialStatus * >( this->giveStatus(gp) );
    bp::object val=ms->giveStateDictionary()[std::to_string(type).c_str()];
    // call parent if we don't have this type in our records
    if(!val) return StructuralMaterial::giveIPValue(answer,gp,type,tStep);
    bp::extract<double> exNum(val);
    bp::extract<FloatArray> exMat(val);
    if(exNum.check()){ answer=FloatArray{exNum()}; return 1; }
    else if(exMat.check()){ answer=exMat(); return 1;} 
    OOFEM_WARNING("Dictionary entry of material state not float or something convertible to FloatArray");
    return 0;
}





void StructuralPythonMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    ///@todo What to do here? Reset dicitonaries? I don't like this function at all.
}


StructuralPythonMaterialStatus :: StructuralPythonMaterialStatus(GaussPoint *gp) :
    StructuralMaterialStatus(gp)
{
}


void StructuralPythonMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    // Copy the temp dict to the equilibrated one
    this->stateDict = this->tempStateDict.copy(); ///@todo Does this suffice? I'm not sure about what happens to references into the dictionary itself. I want a deep copy. / Mikael
}


void StructuralPythonMaterialStatus :: reinitTempStateDictionary()
{
    tempStateDict = stateDict.copy();
}


} // end namespace oofem
