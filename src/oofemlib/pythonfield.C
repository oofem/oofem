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

//#include <Python.h>
#include <pybind11/embed.h>

#include "field.h"
#include "pythonfield.h"
#include "floatarray.h"
#include "timestep.h"

namespace py = pybind11;

namespace oofem {
// REGISTER_Material(PythonField);


/// Constructor.
PythonField :: PythonField(void) : Field(FieldType::FT_Unknown)
{
    //pybind11::print("Hello, World!"); // use the Python API
}

void PythonField :: setFunctionName(std::string functionName){
    this->functionName = functionName;
}


void PythonField :: setModuleName(std::string moduleName){
    this->moduleName = moduleName;
    //this->moduleName.resize(this->moduleName.size()); //remove trailing quotes
}


int PythonField :: evaluateAt(FloatArray &answer, const FloatArray &coords, ValueModeType mode, TimeStep *tStep)
{
    // Do not use the raw CPython API functions Py_Initialize and Py_Finalize as these do not properly handle the lifetime of pybind11’s internal data.
    // https://pybind11.readthedocs.io/en/stable/advanced/embedding.html
    
    // py::scoped_interpreter guard{}; // start the interpreter and keep it alive
    
//     py::initialize_interpreter();
    
    py::module calc = py::module::import(moduleName.c_str());
    py::object result = calc.attr(functionName.c_str())(coords, mode, tStep);
    answer = result.cast<FloatArray>();
//     py::finalize_interpreter();
    
    return 0;
}

int PythonField :: evaluateAt(FloatArray &answer, DofManager *dman, ValueModeType mode, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented");
    return 1;
}



} // end namespace oofem
