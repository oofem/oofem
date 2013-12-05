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

#include "userdefdirichletbc.h"

#include "boundarycondition.h"
#include "timestep.h"
#include "loadtimefunction.h"
#include "verbose.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "dofmanager.h"

#include <Python.h>

namespace oofem {
REGISTER_BoundaryCondition(UserDefDirichletBC);

UserDefDirichletBC :: UserDefDirichletBC(int i, Domain *d) : BoundaryCondition(i, d)
{}


UserDefDirichletBC :: ~UserDefDirichletBC()
{
    Py_DECREF(mpName);
    Py_DECREF(mpModule);
    Py_DECREF(mpFunc);
}


double
UserDefDirichletBC :: give(Dof *dof, ValueModeType mode, TimeStep *stepN)
{
    double factor = this->giveLoadTimeFunction()->evaluate(stepN, mode);
    DofManager *dMan = dof->giveDofManager();


    /*
     * The Python function takes two input arguments:
     *  1) An array with node coordinates
     *  2) The dof id
     */
    int numArgs = 3;

    // Create array with node coordinates
    int dim = dMan->giveCoordinates()->giveSize();
    PyObject *pArgArray = PyList_New(dim);

    PyObject *pArgs = PyTuple_New(numArgs);

    for ( int i = 0; i < dim; i++ ) {
        PyList_SET_ITEM( pArgArray, i, PyFloat_FromDouble( dMan->giveCoordinate(i + 1) ) );
    }

    // PyTuple_SetItem takes over responsibility for objects passed
    // to it -> no DECREF
    PyTuple_SetItem(pArgs, 0, pArgArray);


    // Dof number
    PyObject *pValDofNum = PyLong_FromLong( dof->giveDofID() );
    PyTuple_SetItem(pArgs, 1, pValDofNum);


    // Time
    PyObject *pTargetTime = PyFloat_FromDouble( stepN->giveTargetTime() );
    PyTuple_SetItem(pArgs, 2, pTargetTime);

    // Value returned from the Python function
    PyObject *pRetVal;

    if ( PyCallable_Check(mpFunc) ) {
        pRetVal = PyObject_CallObject(mpFunc, pArgs);
    } else {
        OOFEM_ERROR("UserDefDirichletBC :: give: Python function is not callable.");
    }

    // Get return value
    double retVal = 0.0;
    if ( pRetVal != NULL ) {
        retVal = PyFloat_AsDouble(pRetVal);
    } else   {
        OOFEM_ERROR("UserDefDirichletBC :: give: Failed to fetch Python return value.");
    }

    // Decrement reference count on pointers
    Py_DECREF(pArgs);
    Py_DECREF(pRetVal);

    return retVal * factor;
}


IRResultType
UserDefDirichletBC :: initializeFrom(InputRecord *ir)
// Sets up the dictionary where the receiver stores the conditions it
// imposes.
{
    GeneralBoundaryCondition :: initializeFrom(ir);

    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->mFileName, _IFT_UserDefDirichletBC_filename);

    // Import Python file
    mpName = PyString_FromString( this->mFileName.c_str() );
    mpModule = PyImport_Import(mpName);

    if ( mpModule != NULL ) {
        // Load and call Python function
        mpFunc = PyObject_GetAttrString(mpModule, "giveUserDefBC");
    }

    Py_INCREF(mpName);
    Py_INCREF(mpModule);
    Py_INCREF(mpFunc);

    return IRRT_OK;
}


void
UserDefDirichletBC :: giveInputRecord(DynamicInputRecord &input)
{
    GeneralBoundaryCondition :: giveInputRecord(input);
    input.setField(this->mFileName, _IFT_UserDefDirichletBC_filename);
}


void
UserDefDirichletBC :: setPrescribedValue(double s)
{
    values.zero();
    values.add(s);
}


void
UserDefDirichletBC :: scale(double s)
{
    values.times(s);
}
} // end namespace oofem
