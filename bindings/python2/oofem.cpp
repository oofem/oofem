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
#include "timestep.h"
#include "domain.h"
#include "engngm.h"
#include "element.h"
#include "field.h"
#include "util.h"

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
    
py::class_<oofem::TimeStep>(m, "TimeStep")
    .def("giveNumber", &oofem::TimeStep::giveNumber)
    .def("giveTargetTime", &oofem::TimeStep::giveTargetTime)
    .def("giveIntrinsicTime", &oofem::TimeStep::giveIntrinsicTime)
    .def("giveTimeIncrement", &oofem::TimeStep::giveTimeIncrement)
    .def("setTimeIncrement", &oofem::TimeStep::setTimeIncrement)
    .def("setTime", &oofem::TimeStep::setTime)
    .def("setTargetTime", &oofem::TimeStep::setTargetTime)
    .def("setIntrinsicTime", &oofem::TimeStep::setIntrinsicTime)
    .def("isNotTheLastStep", &oofem::TimeStep::isNotTheLastStep)
    .def("isTheFirstStep", &oofem::TimeStep::isTheFirstStep)
    .def("isTheCurrentTimeStep", &oofem::TimeStep::isTheCurrentTimeStep)
    ;

    py::class_<oofem::EngngModel>(m, "EngngModel")
        .def("giveDomain", &oofem::EngngModel::giveDomain, py::return_value_policy::reference)
        .def("setDomain", &oofem::EngngModel::setDomain, py::keep_alive<1, 3>())
        .def("giveNimberOfDomains", &oofem::EngngModel::giveNumberOfDomains)
        .def("terminateAnalysis", &oofem::EngngModel::terminateAnalysis)
        .def("solveYourself", &oofem::EngngModel::solveYourself)
        .def("solveYourselfAt", &oofem::EngngModel::solveYourselfAt)
        .def("terminate",&oofem::EngngModel::terminate)
        .def("giveField", &oofem::EngngModel::giveField)
        .def("giveCurrentStep", &oofem::EngngModel::giveCurrentStep, py::return_value_policy::reference)
        .def("givePreviousStep", &oofem::EngngModel::givePreviousStep, py::return_value_policy::reference)
        .def("giveNextStep", &oofem::EngngModel::giveNextStep, py::return_value_policy::reference)
        .def("giveNumberOfSteps", &oofem::EngngModel::giveNumberOfSteps)
        //.def("giveUnknownComponent", &oofem::EngngModel::giveUnknownComponent)
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
    //m.def("InstanciateProblem", &oofem::InstanciateProblem);

    py::enum_<oofem::FieldType>(m, "FieldType")
        .value("FT_Unknown", oofem::FieldType::FT_Unknown) 
        .value("FT_Velocity", oofem::FieldType::FT_Velocity)
        .value("FT_Displacements", oofem::FieldType::FT_Displacements)
        .value("FT_VelocityPressure", oofem::FieldType::FT_VelocityPressure)
        .value("FT_Pressure", oofem::FieldType::FT_Pressure)
        .value("FT_Temperature", oofem::FieldType::FT_Temperature)
        .value("FT_HumidityConcentration", oofem::FieldType::FT_HumidityConcentration)
        .value("FT_TransportProblemUnknowns", oofem::FieldType::FT_TransportProblemUnknowns)
        .value("FT_TemperatureAmbient", oofem::FieldType::FT_TemperatureAmbient)
    ;

}