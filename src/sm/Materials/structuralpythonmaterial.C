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
    StructuralMaterial(n, d),
    smallDef(NULL),
    smallDefTangent(NULL),
    largeDef(NULL),
    largeDefTangent(NULL)
{}

StructuralPythonMaterial :: ~StructuralPythonMaterial()
{
    if ( mpModule ) Py_DECREF(mpModule);
}

IRResultType StructuralPythonMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, this->moduleName, _IFT_StructuralPythonMaterial_moduleName);

    // Import Python file
    PyObject *mpName = PyString_FromString( this->moduleName.c_str() );
    this->mpModule = PyImport_Import(mpName);
    Py_DECREF(mpName);

    if ( mpModule != NULL ) {
        // Load and call Python function
        smallDef = PyObject_GetAttrString(mpModule, "computeStress");
        largeDef = PyObject_GetAttrString(mpModule, "computePK1Stress");
        smallDefTangent = PyObject_GetAttrString(mpModule, "computeStressTangent");
        largeDefTangent = PyObject_GetAttrString(mpModule, "computePK1StressTangent");
        if ( smallDefTangent == NULL && smallDef != NULL) {
            OOFEM_WARNING("Using numerical tangent for small deformations");
        }
        if ( largeDefTangent == NULL && largeDef != NULL) {
            OOFEM_WARNING("Using numerical tangent for large deformations");
        }
        if ( smallDef == NULL && largeDef == NULL ) {
            OOFEM_WARNING("No functions for either small or large deformations supplied. Are you sure the functions are named correctly?");
            return IRRT_BAD_FORMAT;
        }
    } else {
        OOFEM_WARNING("mpModule == NULL for module name %s", this->moduleName.c_str());
        return IRRT_BAD_FORMAT;
    }
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
    return new StructuralPythonMaterialStatus(this->giveDomain(), gp);
}

void StructuralPythonMaterial :: callStressFunction(PyObject *func, const FloatArray &oldStrain, const FloatArray &oldStress, const FloatArray &strain, FloatArray &stress, PyObject *stateDict, PyObject *tempStateDict, TimeStep *tStep) const
{
    if ( !PyCallable_Check(func) ) {
        OOFEM_ERROR("Python function is not callable.");
    }

    // Build the argument list;
    double size = strain.giveSize();
    PyObject *pArgs = PyTuple_New(5);

    PyObject *pArgOldStrain = PyList_New(size);
    PyObject *pArgOldStress = PyList_New(size);
    PyObject *pArgStrain = PyList_New(size);

    for ( int i = 0; i < size; i++ ) {
        PyList_SET_ITEM( pArgOldStrain, i, PyFloat_FromDouble( oldStrain[i] ) );
        PyList_SET_ITEM( pArgOldStress, i, PyFloat_FromDouble( oldStress[i] ) );
        PyList_SET_ITEM( pArgStrain, i, PyFloat_FromDouble( strain[i] ) );
    }

    // PyTuple_SetItem takes over responsibility for objects passed to it -> no DECREF
    PyTuple_SetItem(pArgs, 0, pArgOldStrain);
    PyTuple_SetItem(pArgs, 1, pArgOldStress);
    PyTuple_SetItem(pArgs, 2, pArgStrain);
    // Internal state variables
    Py_INCREF(stateDict); ///@todo Verify this; we don't want pArgs to take over ownership, so is this the right thing to do? / Mikael
    PyTuple_SetItem(pArgs, 3, stateDict);
    Py_INCREF(tempStateDict);
    PyTuple_SetItem(pArgs, 4, tempStateDict);
    // Time
    PyTuple_SetItem(pArgs, 5, PyFloat_FromDouble( tStep->giveTargetTime() ));
    // Call the function;
    PyObject *retVal = PyObject_CallObject(func, pArgs);

    // Convert function output back to C++ form:
    stress.resize(size);
    for ( int i = 0; i < size; i++ ) {
        PyObject *val = PyList_GET_ITEM(retVal, i);
        stress[i] = PyFloat_AS_DOUBLE(val);
        //Py_DECREF(val);
    }

    Py_DECREF(retVal);
    Py_DECREF(pArgs);
}

void StructuralPythonMaterial :: callTangentFunction(FloatMatrix &answer, PyObject *func, const FloatArray &strain, const FloatArray &stress, PyObject *stateDict, PyObject *tempStateDict, TimeStep *tStep) const
{
    if ( !PyCallable_Check(func) ) {
        OOFEM_ERROR("Python function is not callable.");
    }

    // Build the argument list;
    double size = strain.giveSize();
    PyObject *pArgs = PyTuple_New(4);

    PyObject *pArgStrain = PyList_New(size);
    PyObject *pArgStress = PyList_New(size);

    for ( int i = 0; i < size; i++ ) {
        PyList_SET_ITEM( pArgStrain, i, PyFloat_FromDouble( strain[i] ) );
        PyList_SET_ITEM( pArgStress, i, PyFloat_FromDouble( stress[i] ) );
    }
    // PyTuple_SetItem takes over responsibility for objects passed to it -> no DECREF
    PyTuple_SetItem(pArgs, 0, pArgStrain);
    PyTuple_SetItem(pArgs, 1, pArgStress);
    // Internal state variables
    Py_INCREF(stateDict); ///@todo Verify this; we don't want pArgs to take over ownership, so is this the right thing to do? / Mikael
    PyTuple_SetItem(pArgs, 2, stateDict);
    Py_INCREF(tempStateDict); ///@todo Verify this; we don't want pArgs to take over ownership, so is this the right thing to do? / Mikael
    PyTuple_SetItem(pArgs, 3, tempStateDict);
    // Time
    PyTuple_SetItem(pArgs, 4, PyFloat_FromDouble( tStep->giveTargetTime() ));
    // Call the function;
    PyObject *retVal = PyObject_CallObject(func, pArgs);

    if ( retVal == NULL ) {
        OOFEM_ERROR("bad return value (null) from PyObject_CallObject");
    }
    // Convert function output back to C++ form:
    answer.resize(size, size);
    for ( int i = 0; i < size; i++ ) {
        PyObject *row = PyList_GET_ITEM(retVal, i); // Get item just borrows the reference, don't DECREF
        for ( int j = 0; j < size; j++ ) {
            PyObject *val = PyList_GET_ITEM(row, j);
            answer(i, j) = PyFloat_AS_DOUBLE(val);
        }
    }

    Py_DECREF(retVal);
    Py_DECREF(pArgs);
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
    std :: string s = std :: to_string( type );
    PyObject *val = PyDict_GetItemString(ms->giveStateDictionary(), s.c_str());
    if ( val ) {
        
        if ( PyFloat_Check(val) != 0 ) {
            answer = FloatArray{PyFloat_AS_DOUBLE(val)};
        } else if ( PyList_Check(val) != 0 ) {
            // Convert function output back to C++ form:
            int size = PyList_Size(val);
            answer.resize(size);
            for ( int i = 0; i < size; i++ ) {
                PyObject *val = PyList_GET_ITEM(val, i);
                answer[i] = PyFloat_AS_DOUBLE(val);
                //Py_DECREF(val);
            }
            return 1;
        } else {
            OOFEM_WARNING("Dictionary entry of material state is not a list (only lists supported for now)");
            return 0;
        }
    }
    return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
}





void StructuralPythonMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    ///@todo What to do here? Reset dicitonaries? I don't like this function at all.
}


StructuralPythonMaterialStatus :: StructuralPythonMaterialStatus(Domain *d, GaussPoint *gp) :
    StructuralMaterialStatus(0, d, gp)
{
    this->stateDict = PyDict_New();
    this->tempStateDict = PyDict_New();
}

StructuralPythonMaterialStatus :: ~StructuralPythonMaterialStatus()
{
    Py_DECREF(this->stateDict);
    Py_DECREF(this->tempStateDict);
}

void StructuralPythonMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    // Copy the temp dict to the equilibrated one
    auto oldDict = this->stateDict;
    this->stateDict = PyDict_Copy(this->tempStateDict); ///@todo Does this suffice? I'm not sure about what happens to references into the dictionary itself. I want a deep copy. / Mikael
    Py_DECREF(oldDict);
}


void StructuralPythonMaterialStatus :: reinitTempStateDictionary()
{
    Py_DECREF(this->tempStateDict);
    this->tempStateDict = PyDict_Copy(this->stateDict);
}


} // end namespace oofem
