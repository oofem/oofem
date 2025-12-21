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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "pythonexpression.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "error.h"

#include <sstream>
#include <iostream>
#include <fstream>

// Defines the name for the return variable;
#define RETURN_VARIABLE "ret"

namespace oofem {
REGISTER_Function(PythonExpression);

PythonExpression :: PythonExpression(int n, Domain *d) : Function(n, d) { }

PythonExpression :: ~PythonExpression()
{
    ///@todo Is this right?
    //Py_DECREF(this->main_dict);
}

void
PythonExpression :: initializeFrom(InputRecord &ir)
{
    Function :: initializeFrom(ir);

    // Check if the f expression is given
    if (ir.hasField(_IFT_PythonExpression_f)) {
        IR_GIVE_FIELD(ir, this->fExpression, _IFT_PythonExpression_f);
    } else {
        std::string path;
        IR_GIVE_FIELD(ir, path, _IFT_PythonExpression_ffile);
        this->readFile2String(path, this->fExpression);
    }
    this->f = Py_CompileString(fExpression.c_str(), "<internal_f>", Py_file_input);
    if ( this->f == nullptr) {
       PyErr_Print();
    }

    // Check if the dfdt expression is given
    if (ir.hasField(_IFT_PythonExpression_dfdt)) {
        IR_GIVE_FIELD(ir, this->dfdtExpression, _IFT_PythonExpression_dfdt);
    } else if (ir.hasField(_IFT_PythonExpression_dfdtfile)) {
        std::string path;
        IR_GIVE_OPTIONAL_FIELD(ir, path, _IFT_PythonExpression_dfdtfile);
        this->readFile2String(path, this->dfdtExpression);
    }
    this->dfdt = Py_CompileString(dfdtExpression.c_str(), "<internal_dfdt>", Py_file_input);
    if ( this->dfdt == nullptr) {
       PyErr_Print();
    }

    // Check if the d2fdt2 expression is given
    if (ir.hasField(_IFT_PythonExpression_d2fdt2)) {
        IR_GIVE_FIELD(ir, this->d2fdt2Expression, _IFT_PythonExpression_d2fdt2);
    } else if (ir.hasField(_IFT_PythonExpression_d2fdt2file)) {
        std::string path;
        IR_GIVE_OPTIONAL_FIELD(ir, path, _IFT_PythonExpression_d2fdt2file);
        this->readFile2String(path, this->d2fdt2Expression);
    }
    this->d2fdt2 = Py_CompileString(d2fdt2Expression.c_str(), "<internal_d2fdt2>", Py_file_input);
    if ( this->d2fdt2 == nullptr) {
       PyErr_Print();
    }

    ///@todo Check this stuff; Is this OK to do? We need a way to fetch the global dictionary..
    if ( !main_dict ) {
        PyObject *main_module = PyImport_ImportModule("__main__");
        if (main_module != NULL) {
            this->main_dict = PyModule_GetDict(main_module);
            Py_DECREF(main_module);  // Decrease reference count
        } else {
            // Handle error: could not import __main__ module
            OOFEM_WARNING("Could not import __main__ module");
            this->main_dict = PyDict_New();
        }
    }
}


void
PythonExpression :: giveInputRecord(DynamicInputRecord &input)
{
    Function :: giveInputRecord(input);
    input.setField(this->fExpression, _IFT_PythonExpression_f);
    input.setField(this->dfdtExpression, _IFT_PythonExpression_dfdt);
    input.setField(this->d2fdt2Expression, _IFT_PythonExpression_d2fdt2);
}


PyObject *
PythonExpression :: getDict(const std :: map< std :: string, FunctionArgument > &valDict)
{
    PyObject *local_dict = PyDict_New();
    for ( const auto &named_arg: valDict ) {
        const FunctionArgument &arg = named_arg.second;
        PyObject *tmp;
        if ( arg.type == FunctionArgument :: FAT_double ) {
            tmp = PyFloat_FromDouble(arg.val0);
        } else if ( arg.type == FunctionArgument :: FAT_FloatArray ) {
            tmp = PyList_New( arg.val1.giveSize() );
            for ( int i = 0; i < arg.val1.giveSize(); ++i ) {
                PyList_SET_ITEM( tmp, i, PyFloat_FromDouble( arg.val1[i] ) );
            }
        } else if ( arg.type == FunctionArgument :: FAT_int ) {
            tmp = PyLong_FromLong(arg.val2);
        } else if ( arg.type == FunctionArgument :: FAT_IntArray ) {
            tmp = PyList_New( arg.val3.giveSize() );
            for ( int i = 0; i < arg.val3.giveSize(); ++i ) {
                PyList_SET_ITEM( tmp, i, PyLong_FromLong( arg.val3[i] ) );
            }
        } else {
            tmp = NULL;
            OOFEM_ERROR("Unsupported FunctionArgumentType");
        }
        PyDict_SetItemString(local_dict, named_arg.first.c_str(), tmp);
    }
    return local_dict;
}


void
PythonExpression :: getArray(FloatArray &answer, PyObject **func, const std :: map< std :: string, FunctionArgument > &valDict)
{
#ifdef _OPENMP
#pragma omp critical
#endif
{
    PyObject *local_dict = getDict(valDict);
    PyObject *dummy = PyEval_EvalCode( *func, main_dict, local_dict );

    if ( dummy == nullptr ) {
        PyErr_Print();
    }
    PyObject *ret = PyDict_GetItemString(local_dict, RETURN_VARIABLE);
    if ( PyList_Check(ret) ) {
        int size = PyList_GET_SIZE(ret);
        answer.resize(size);
        for ( int i = 0; i < size; ++i ) {
            answer(i) = this->pyObj2double( PyList_GET_ITEM(ret, i) );
        }
    } else {
        answer = {this->pyObj2double(ret) };
    }

    Py_DECREF(local_dict);
    Py_DECREF(dummy);
    Py_DECREF(ret);
}
}


void
PythonExpression :: evaluate(FloatArray &answer, const std :: map< std :: string, FunctionArgument > &valDict, GaussPoint *gp, double param)
{
    this->getArray(answer, &this->f, valDict);
}


void
PythonExpression :: evaluateVelocity(FloatArray &answer, const std :: map< std :: string, FunctionArgument > &valDict)
{
    this->getArray(answer, &this->dfdt, valDict);
}


void
PythonExpression :: evaluateAcceleration(FloatArray &answer, const std :: map< std :: string, FunctionArgument > &valDict)
{
    this->getArray(answer, &this->d2fdt2, valDict);
}


double
PythonExpression :: getScalar(PyObject *func, double time)
{
    PyObject *local_dict = PyDict_New();
    PyDict_SetItemString( local_dict, "t", PyFloat_FromDouble(time) );
    PyObject *dummy = PyEval_EvalCode( func, main_dict, local_dict );
    double val = 0.;
    PyObject *ret = PyDict_GetItemString(local_dict, RETURN_VARIABLE);
    if ( PyNumber_Check(ret) ) {
        val = pyObj2double(ret);
    } else if ( PyList_Check(ret) ) {
        if ( PyList_GET_SIZE(ret) != 1 ) {
            OOFEM_ERROR("Result from python is not a real float!");
        } else {
            val = pyObj2double( PyList_GET_ITEM(ret, 0) );
        }
    } else {
        OOFEM_ERROR("Result from python is not a real float!");
    }

    Py_DECREF(local_dict);
    Py_DECREF(dummy);
    Py_DECREF(ret);
    return val;
}


double PythonExpression :: evaluateAtTime(double time)
{
    return this->getScalar(this->f, time);
}

double PythonExpression :: evaluateVelocityAtTime(double time)
{
    return this->getScalar(this->dfdt, time);
}


double PythonExpression :: evaluateAccelerationAtTime(double time)
{
    return this->getScalar(this->d2fdt2, time);
}


double PythonExpression::pyObj2double(PyObject *obj) {
    if (PyNumber_Check(obj)) {
        PyObject *float_obj = PyNumber_Float(obj);
        if (float_obj != NULL) {
            double result = PyFloat_AsDouble(float_obj);
            Py_DECREF(float_obj);  // Decrease reference count
            return result;
        }
    }
    // Handle error: object is not a number or conversion failed
    return -1.0;  // Or some other error indicator
}


void PythonExpression::readFile2String(const std::string &path, std::string &content) {
    std::ifstream file(path);
    if (file.is_open()) {
        std::stringstream buffer;
        buffer << file.rdbuf();
        content = buffer.str();
        file.close();
    } else {
        OOFEM_ERROR("Could not open file %s", path.c_str());
    }
}

} // end namespace oofem
