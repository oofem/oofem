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
#include "floatmatrix.h"
#include "classfactory.h"
#include "domain.h"
#include "engngm.h"
#include "element.h"


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

     py::class_<oofem::FloatMatrix>(m, "FloatMatrix")
        .def(py::init<>())
        .def(py::init<int,int>())
        .def("printYourself", (void (oofem::FloatMatrix::*)() const) &oofem::FloatMatrix::printYourself, "Prints receiver")
        .def("printYourself", (void (oofem::FloatMatrix::*)(const std::string &) const) &oofem::FloatMatrix::printYourself, "Prints receiver")
        .def("pY", &oofem::FloatMatrix::pY)
        .def("__setitem__", [](oofem::FloatMatrix &s, py::tuple indx, double v) {
            s(((py::handle)(indx[0])).cast<int>(), ((py::handle)(indx[1])).cast<int>()) = v;
        })
        .def("__getitem__", [](const oofem::FloatMatrix &s, py::tuple indx) {
            if (((py::handle)(indx[0])).cast<int>() >= s.giveNumberOfRows()) throw py::index_error();
            if (((py::handle)(indx[1])).cast<int>() >= s.giveNumberOfColumns()) throw py::index_error();
            return s(((py::handle)(indx[0])).cast<int>(), ((py::handle)(indx[1])).cast<int>());
        })
        ;
    
    py::class_<oofem::EngngModel>(m, "EngngModel")
        .def("giveDomain", &oofem::EngngModel::giveDomain, py::return_value_policy::reference)
        ;

    py::class_<oofem::Domain>(m, "Domain")
    ;

    py::class_<oofem::Element>(m, "Element")
    ;

    py::class_<oofem::ClassFactory>(m, "ClassFactory")
        .def("createElement", &oofem::ClassFactory::createElement)
        .def("createEngngModel", &oofem::ClassFactory::createEngngModel)
        ;

    m.def("getClassFactory", &oofem::GiveClassFactory, py::return_value_policy::reference);

}