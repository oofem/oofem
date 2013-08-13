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

#ifdef BOOST_PYTHON
#include <Python.h>
#endif

namespace oofem {


REGISTER_BoundaryCondition( UserDefDirichletBC );

#ifdef BOOST_PYTHON
PythonInitializer mPythonInitializer;
#endif

double UserDefDirichletBC :: give(Dof *dof, ValueModeType mode, TimeStep *stepN)
{
#ifdef BOOST_PYTHON
    double factor = this->giveLoadTimeFunction()->evaluate(stepN, mode);
    DofManager *dMan = dof->giveDofManager();


    // Import Python file
    PyObject *pName = PyString_FromString( mFileName.c_str() );
    PyObject *pModule = PyImport_Import( pName );


    /*
     * Create input arguments to the Python function:
     * (x, y, z, dofNum)
     */
    int numArgs = 4;
    PyObject *pArgs = PyTuple_New(numArgs);

    // X-coordinate
    PyObject *pValX = PyFloat_FromDouble(dMan->giveCoordinate(1));
    PyTuple_SetItem(pArgs, 0, pValX);

    // Y-coordinate
    PyObject *pValY = PyFloat_FromDouble(dMan->giveCoordinate(2));
    PyTuple_SetItem(pArgs, 1, pValY);

    // Z-coordinate
    PyObject *pValZ = PyFloat_FromDouble(dMan->giveCoordinate(3));
    PyTuple_SetItem(pArgs, 2, pValZ);


    // Dof number
    PyObject *pValDofNum = PyLong_FromLong(dof->giveNumber());
    PyTuple_SetItem(pArgs, 3, pValDofNum);



    PyObject *pRetVal;

    if (pModule != NULL) {
    	// Load and call Python function
    	PyObject *pFunc = PyObject_GetAttrString(pModule, "giveUserDefBC");

    	if(pFunc != NULL) {
    	    if (PyCallable_Check(pFunc)) {
    	    	pRetVal = PyObject_CallObject(pFunc, pArgs);
    	    } else {
    	    	OOFEM_ERROR("UserDefDirichletBC :: give: Python function is not callable.");
    	    }
    	}
    	else {
	    	OOFEM_ERROR("UserDefDirichletBC :: give: Failed to load Python function.");
    	}
    }
    else {
    	OOFEM_ERROR("UserDefDirichletBC :: give: Failed to load Python module.");
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
    Py_DECREF(pName);
    Py_DECREF(pModule);
    Py_DECREF(pArgs);
    Py_DECREF(pRetVal);


    return retVal*factor;
#else
	OOFEM_ERROR("UserDefDirichletBC :: give requires BOOST_PYTHON.");
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

#ifdef BOOST_PYTHON
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
