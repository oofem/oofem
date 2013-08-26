/*
 * userdefdirichletbc.C
 *
 *  Created on: Aug 7, 2013
 *      Author: svennine
 */

#include "userdefdirichletbc.h"

#include "boundarycondition.h"
#include "timestep.h"
#include "loadtimefunction.h"
#include "verbose.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"

#include "dofmanager.h"

#ifdef USERDEFDIRICHLETBC
#include <Python.h>
#endif

namespace oofem {


REGISTER_BoundaryCondition( UserDefDirichletBC );

#ifdef USERDEFDIRICHLETBC
PythonInitializer mPythonInitializer;
#endif

UserDefDirichletBC :: UserDefDirichletBC(int i, Domain *d) : BoundaryCondition(i, d)
{

}

/// Destructor
UserDefDirichletBC :: ~UserDefDirichletBC()
{
#ifdef USERDEFDIRICHLETBC
    Py_DECREF(mpName);
    Py_DECREF(mpModule);
    Py_DECREF(mpFunc);
#endif
}

double UserDefDirichletBC :: give(Dof *dof, ValueModeType mode, TimeStep *stepN)
{
#ifdef USERDEFDIRICHLETBC
    double factor = this->giveLoadTimeFunction()->evaluate(stepN, mode);
    DofManager *dMan = dof->giveDofManager();


    /**
     * The Python function takes two input arguments:
     * 	1) An array with node coordinates
     * 	2) The dof id
     */
    int numArgs = 2;

    // Create array with node coordinates
    int dim = dMan->giveCoordinates()->giveSize();
    PyObject *pArgArray = PyList_New(dim);

    PyObject *pArgs = PyTuple_New(numArgs);

    for (int i = 0; i < dim; i++) {
        PyList_SET_ITEM(pArgArray, i, PyFloat_FromDouble( dMan->giveCoordinate(i+1) ));
    }
    // PyTuple_SetItem takes over responsibility for objects passed
    // to it -> no DECREF
    PyTuple_SetItem(pArgs, 0, pArgArray);


    // Dof number
    PyObject *pValDofNum = PyLong_FromLong(dof->giveDofID());
    PyTuple_SetItem(pArgs, 1, pValDofNum);


    // Value returned from the Python function
    PyObject *pRetVal;

    if (PyCallable_Check(mpFunc)) {
    	pRetVal = PyObject_CallObject(mpFunc, pArgs);
    } else {
    	OOFEM_ERROR("UserDefDirichletBC :: give: Python function is not callable.");
    }

    // Get return value
    double retVal = 0.0;
    if(pRetVal != NULL) {
    	retVal = PyFloat_AsDouble(pRetVal);
    }
    else {
    	OOFEM_ERROR("UserDefDirichletBC :: give: Failed to fetch Python return value.");
    }

    // Decrement reference count on pointers
    Py_DECREF(pArgs);
    Py_DECREF(pRetVal);


    return retVal*factor;
#else
	OOFEM_ERROR("UserDefDirichletBC :: give requires USERDEFDIRICHLETBC to be defined (to enable Python).");
    return 0.0;
#endif
}


IRResultType
UserDefDirichletBC :: initializeFrom(InputRecord *ir)
// Sets up the dictionary where the receiver stores the conditions it
// imposes.
{
    GeneralBoundaryCondition :: initializeFrom(ir);

    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, mFileName, _IFT_UserDefDirichletBC_filename);

#ifdef USERDEFDIRICHLETBC
    // Import Python file
    mpName = PyString_FromString( mFileName.c_str() );
    mpModule = PyImport_Import( mpName );

    if (mpModule != NULL) {
    	// Load and call Python function
    	mpFunc = PyObject_GetAttrString(mpModule, "giveUserDefBC");
    }

    Py_INCREF(mpName);
    Py_INCREF(mpModule);
    Py_INCREF(mpFunc);

#endif

    return IRRT_OK;
}


void
UserDefDirichletBC :: giveInputRecord(DynamicInputRecord &input)
{
//    GeneralBoundaryCondition :: giveInputRecord(input);
//    input.setField(this->values, _IFT_BoundaryCondition_values);
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

#ifdef USERDEFDIRICHLETBC
PythonInitializer :: PythonInitializer()
{
    Py_Initialize();


    // Adding . to the system path allows us to run Python
    // functions stored in the working directory.
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");
}

PythonInitializer :: ~PythonInitializer()
{
    Py_Finalize();
}
#endif

} // end namespace oofem
