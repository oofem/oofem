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

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "floatarray.h"


PYBIND11_MODULE(oofempy, m) {
    m.doc() = "oofem python bindings module"; // optional module docstring

    
    py::class_<oofem::FloatArray>(m, "FloatArray")
        .def(py::init<int>(), py::arg("n")=0)
        .def("printYourself", (void (oofem::FloatArray::*)() const) &oofem::FloatArray::printYourself, "Prints receiver")
        .def("printYourself", (void (oofem::FloatArray::*)(const std::string &) const) &oofem::FloatArray::printYourself, "Prints receiver")
        .def("pY", &oofem::FloatArray::pY)
        .def("__setitem__", [](oofem::FloatArray &s, size_t i, double v) {
            s[i] = v;
        })
        .def("__getitem__", [](const oofem::FloatArray &s, size_t i) {
            if (i >= (size_t) s.giveSize()) throw py::index_error();
            return s[i];
        })
        ;
      

}