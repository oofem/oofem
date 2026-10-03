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

namespace {
/**
 * RAII acquisition of the global interpreter lock.
 *
 * PyGILState_Ensure() succeeds from any thread and nests safely, so the same guard serves
 * the main thread during instanciation and the worker threads of a parallel assembly loop.
 * Declare it before any PyRef, so that the references are still released under the lock.
 */
class GILGuard
{
public:
    GILGuard()
    {
        if ( !Py_IsInitialized() ) {
            OOFEM_ERROR("The Python interpreter has not been initialized");
        }
        state = PyGILState_Ensure();
    }
    ~GILGuard() { PyGILState_Release(state); }
    GILGuard(const GILGuard &) = delete;
    GILGuard &operator=(const GILGuard &) = delete;

private:
    PyGILState_STATE state;
};

/// Owning handle for a Python reference, released on scope exit, including while unwinding.
class PyRef
{
public:
    explicit PyRef(PyObject * o = nullptr) : obj(o) { }
    ~PyRef() { Py_XDECREF(obj); }
    PyRef(const PyRef &) = delete;
    PyRef &operator=(const PyRef &) = delete;
    operator PyObject *() const { return obj; }
    /// Hands over ownership to the caller.
    PyObject *release() { PyObject *o = obj; obj = nullptr; return o; }

private:
    PyObject *obj;
};
} // end anonymous namespace

/// Compiles an expression, reporting a syntax error rather than leaving a null code object behind.
static PyObject *
compileExpression(const std :: string &expression, const char *name)
{
    std :: string filename = std :: string("<internal_") + name + ">";
    PyObject *code = Py_CompileString( expression.c_str(), filename.c_str(), Py_file_input );
    if ( code == nullptr ) {
        PyErr_Print();
        OOFEM_ERROR( "Could not compile the \"%s\" expression \"%s\"", name, expression.c_str() );
    }
    return code;
}

PythonExpression :: PythonExpression(int n, Domain *d) : Function(n, d) { }

PythonExpression :: ~PythonExpression()
{
    // Functions are owned by the domain and may outlive Py_Finalize(), after which
    // touching a Python object is undefined; the references are deliberately left
    // for the interpreter teardown to reclaim.
}

void
PythonExpression :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    Function :: initializeFrom(ir);

    GILGuard gil;

    // Check if the f expression is given
    if (ir->hasField(_IFT_PythonExpression_f)) {
        IR_GIVE_FIELD(ir, this->fExpression, _IFT_PythonExpression_f);
    } else {
        std::string path;
        IR_GIVE_FIELD(ir, path, _IFT_PythonExpression_ffile);
        this->readFile2String(path, this->fExpression);
    }
    this->f = compileExpression(this->fExpression, _IFT_PythonExpression_f);

    // Check if the dfdt expression is given
    if (ir->hasField(_IFT_PythonExpression_dfdt)) {
        IR_GIVE_FIELD(ir, this->dfdtExpression, _IFT_PythonExpression_dfdt);
    } else if (ir->hasField(_IFT_PythonExpression_dfdtfile)) {
        std::string path;
        IR_GIVE_OPTIONAL_FIELD(ir, path, _IFT_PythonExpression_dfdtfile);
        this->readFile2String(path, this->dfdtExpression);
    }
    this->dfdt = compileExpression(this->dfdtExpression, _IFT_PythonExpression_dfdt);

    // Check if the d2fdt2 expression is given
    if (ir->hasField(_IFT_PythonExpression_d2fdt2)) {
        IR_GIVE_FIELD(ir, this->d2fdt2Expression, _IFT_PythonExpression_d2fdt2);
    } else if (ir->hasField(_IFT_PythonExpression_d2fdt2file)) {
        std::string path;
        IR_GIVE_OPTIONAL_FIELD(ir, path, _IFT_PythonExpression_d2fdt2file);
        this->readFile2String(path, this->d2fdt2Expression);
    }
    this->d2fdt2 = compileExpression(this->d2fdt2Expression, _IFT_PythonExpression_d2fdt2);

    ///@todo Check this stuff; Is this OK to do? We need a way to fetch the global dictionary..
    if ( !main_dict ) {
        PyObject *main_module = PyImport_ImportModule("__main__");
        if (main_module != NULL) {
            // PyModule_GetDict returns a borrowed reference; take our own, as the
            // dictionary is kept for the lifetime of this object.
            this->main_dict = PyModule_GetDict(main_module);
            Py_XINCREF(this->main_dict);
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
    PyRef local_dict( PyDict_New() );
    for ( const auto &named_arg: valDict ) {
        const FunctionArgument &arg = named_arg.second;
        PyObject *value = nullptr;
        if ( arg.type == FunctionArgument :: FAT_double ) {
            value = PyFloat_FromDouble(arg.val0);
        } else if ( arg.type == FunctionArgument :: FAT_FloatArray ) {
            value = PyList_New( arg.val1.giveSize() );
            for ( int i = 0; i < arg.val1.giveSize(); ++i ) {
                PyList_SET_ITEM( value, i, PyFloat_FromDouble( arg.val1[i] ) );
            }
        } else if ( arg.type == FunctionArgument :: FAT_int ) {
            value = PyLong_FromLong(arg.val2);
        } else if ( arg.type == FunctionArgument :: FAT_IntArray ) {
            value = PyList_New( arg.val3.giveSize() );
            for ( int i = 0; i < arg.val3.giveSize(); ++i ) {
                PyList_SET_ITEM( value, i, PyLong_FromLong( arg.val3[i] ) );
            }
        } else {
            OOFEM_ERROR( "Unsupported FunctionArgumentType for argument \"%s\"", named_arg.first.c_str() );
        }
        // The dictionary takes its own reference; tmp releases ours.
        PyRef tmp(value);
        PyDict_SetItemString( local_dict, named_arg.first.c_str(), tmp );
    }
    return local_dict.release();
}


void
PythonExpression :: getArray(FloatArray &answer, PyObject **func, const std :: map< std :: string, FunctionArgument > &valDict)
{
    // The lock serialises the evaluation; each call builds its own local dictionary, so
    // concurrent evaluations of the same function do not share any state beyond the globals.
    GILGuard gil;

    PyRef local_dict( this->getDict(valDict) );
    PyRef dummy( PyEval_EvalCode( *func, main_dict, local_dict ) );
    if ( !dummy ) {
        PyErr_Print();
        OOFEM_ERROR("Evaluation of the Python expression failed");
    }

    // Borrowed reference owned by local_dict; it must not be released here.
    PyObject *ret = PyDict_GetItemString(local_dict, RETURN_VARIABLE);
    if ( ret == nullptr ) {
        OOFEM_ERROR("Python expression did not assign the return variable \"" RETURN_VARIABLE "\"");
    }

    if ( PyList_Check(ret) ) {
        int size = PyList_GET_SIZE(ret);
        answer.resize(size);
        for ( int i = 0; i < size; ++i ) {
            answer(i) = this->pyObj2double( PyList_GET_ITEM(ret, i) );
        }
    } else {
        answer = {this->pyObj2double(ret) };
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
    GILGuard gil;

    PyRef local_dict( PyDict_New() );
    PyRef t( PyFloat_FromDouble(time) );
    PyDict_SetItemString(local_dict, "t", t); // the dictionary takes its own reference

    PyRef dummy( PyEval_EvalCode( func, main_dict, local_dict ) );
    if ( !dummy ) {
        PyErr_Print();
        OOFEM_ERROR("Evaluation of the Python expression failed");
    }

    // Borrowed reference owned by local_dict; it must not be released here.
    PyObject *ret = PyDict_GetItemString(local_dict, RETURN_VARIABLE);
    if ( ret == nullptr ) {
        OOFEM_ERROR("Python expression did not assign the return variable \"" RETURN_VARIABLE "\"");
    }

    double val;
    if ( PyNumber_Check(ret) ) {
        val = pyObj2double(ret);
    } else if ( PyList_Check(ret) && PyList_GET_SIZE(ret) == 1 ) {
        val = pyObj2double( PyList_GET_ITEM(ret, 0) );
    } else {
        OOFEM_ERROR("Result from python is not a real float!");
    }
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
