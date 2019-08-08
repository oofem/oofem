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
#include <pybind11/operators.h>
namespace py = pybind11;

#include <string>

#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"

#include "femcmpnn.h"
#include "timestep.h"
#include "domain.h"
#include "engngm.h"
#include "dof.h"
#include "dofmanager.h"
#include "element.h"
#include "generalboundarycondition.h"
#include "initialcondition.h"
#include "function.h"
#include "material.h"
#include "crosssection.h"
#include "field.h"
#include "util.h"
#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "valuemodetype.h"
#include "dofiditem.h"
#include "classfactory.h"


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
    
    py::class_<oofem::IntArray>(m, "IntArray")
        .def(py::init<int>(), py::arg("n")=0)
        .def(py::init<const oofem::IntArray&>())
        .def("resize", &oofem::IntArray::resize)
        .def("clear", &oofem::IntArray::clear)
        .def("preallocate", &oofem::IntArray::preallocate)
        .def("followedBy", (void (oofem::IntArray::*)(int , int)) &oofem::IntArray::followedBy, "Appends number to receiver")
        .def("followedBy", (void (oofem::IntArray::*)(const oofem::IntArray& , int)) &oofem::IntArray::followedBy, "Appends array to receiver")
        .def("isEmpty", &oofem::IntArray::isEmpty)
        .def("containsOnlyZeroes", &oofem::IntArray::containsOnlyZeroes)
        .def("containsSorted", &oofem::IntArray::containsSorted)
        .def("insertSorted", &oofem::IntArray::insertSorted)
        .def("eraseSorted", &oofem::IntArray::eraseSorted)
        .def("findFirstIndexOf", &oofem::IntArray::findFirstIndexOf, "Finds index of first occurrence of given value in array")
        .def("contains", &oofem::IntArray::contains)
        .def("sort", &oofem::IntArray::sort)
        .def("zero", &oofem::IntArray::zero)
        .def("pY", &oofem::IntArray::pY)
        .def("__repr__",
            [](const oofem::IntArray &s) {
                std::string a = "<oofempy.IntArray: {";
                for ( int i = 0; i < s.giveSize(); ++i ) {
                    if ( i > 40 ) {
                        a.append("...");
                        break;
                    } else {
                        a.append(std::to_string(s[i]));
                        a.append(", ");
                    }
                }
                a.append("}>");
                return a;
            })
        .def("__setitem__", [](oofem::IntArray &s, size_t i, int v) {
            s[i] = v;
        })
        .def("__getitem__", [](const oofem::IntArray &s, size_t i) {
            if (i >= (size_t) s.giveSize()) throw py::index_error();
            return s[i];
        })
    ;

    py::class_<oofem::DataReader>(m, "DataReader")
    ;

    py::class_<oofem::OOFEMTXTDataReader, oofem::DataReader>(m, "OOFEMTXTDataReader")
        .def(py::init<std::string>())
    ;

    py::class_<oofem::FEMComponent>(m, "FEMComponent")
        .def("giveClassName", &oofem::FEMComponent::giveClassName)
        .def("giveInputRecordName", &oofem::FEMComponent::giveInputRecordName)
        .def("giveDomain", &oofem::FEMComponent::giveDomain, py::return_value_policy::reference)
        .def("setDomain", &oofem::FEMComponent::setDomain)
        .def("giveNumber", &oofem::FEMComponent::giveNumber)
        .def("setNumber", &oofem::FEMComponent::setNumber)
        .def("checkConsistency", &oofem::FEMComponent::checkConsistency)
        .def("printYourself", &oofem::FEMComponent::printYourself)
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
        .def("giveUnknownComponent", &oofem::EngngModel::giveUnknownComponent)
        ;

    py::class_<oofem::Domain>(m, "Domain")
        .def(py::init<int, int, oofem::EngngModel*>())
        .def("giveNumber", &oofem::Domain::giveNumber)
        .def("setNumber", &oofem::Domain::setNumber)
        .def("giveElement", &oofem::Domain::giveElement, py::return_value_policy::reference)
        .def("giveBc", &oofem::Domain::giveBc, py::return_value_policy::reference)
        .def("giveIc", &oofem::Domain::giveIc, py::return_value_policy::reference)
        .def("giveFunction", &oofem::Domain::giveFunction, py::return_value_policy::reference)
        .def("giveMaterial", &oofem::Domain::giveMaterial, py::return_value_policy::reference)
        .def("giveCrossSection", &oofem::Domain::giveCrossSection, py::return_value_policy::reference)
    ;

    py::class_<oofem::Dof>(m, "Dof")
        .def("giveDofManNumber", &oofem::Dof::giveDofManNumber)
        .def("giveDofManager", &oofem::Dof::giveDofManager, py::return_value_policy::reference)
        .def("giveDofManGlobalNumber", &oofem::Dof::giveDofManGlobalNumber)
        .def("giveUnknown", (double (oofem::Dof::*)(oofem::ValueModeType, oofem::TimeStep*)) &oofem::Dof::giveUnknown)
        .def("giveDofID", &oofem::Dof::giveDofID)
        .def("printYourself", &oofem::Dof::printYourself)
    ;

    py::class_<oofem::DofManager, oofem::FEMComponent>(m, "DofManager")
        .def("giveDofWithID", &oofem::DofManager::giveDofWithID, py::return_value_policy::reference)
        .def("giveNumberOfDofs", &oofem::DofManager::giveNumberOfDofs)
        .def("giveUnknownVector", (void (oofem::DofManager::*)(oofem::FloatArray&, const oofem::IntArray&, oofem::ValueModeType, oofem::TimeStep*, bool)) &oofem::DofManager::giveUnknownVector)
        .def("givePrescribedUnknownVector", &oofem::DofManager::givePrescribedUnknownVector)
        .def("giveCoordinates", &oofem::DofManager::giveCoordinates, py::return_value_policy::reference)
        .def("appendDof", &oofem::DofManager::appendDof, py::keep_alive<1, 2>())
        .def("removeDof", &oofem::DofManager::removeDof)
        .def("hasDofID", &oofem::DofManager::hasDofID)
    ;

    py::class_<oofem::Element, oofem::FEMComponent>(m, "Element")
    ;

    py::class_<oofem::GeneralBoundaryCondition, oofem::FEMComponent>(m, "GeneralBoundaryCondition")
    ;

    py::class_<oofem::InitialCondition, oofem::FEMComponent>(m, "InitialCondition")
    ;

    py::class_<oofem::Function, oofem::FEMComponent>(m, "Function")
    ;

    py::class_<oofem::Material, oofem::FEMComponent>(m, "Material")
    ;

    py::class_<oofem::CrossSection, oofem::FEMComponent>(m, "CrossSection")
    ;


  

    py::class_<oofem::ClassFactory>(m, "ClassFactory")
        .def("createElement", &oofem::ClassFactory::createElement)
        .def("createEngngModel", &oofem::ClassFactory::createEngngModel)
        ;

    m.def("getClassFactory", &oofem::GiveClassFactory, py::return_value_policy::reference);
    m.def("InstanciateProblem", &oofem::InstanciateProblem);
    //std::unique_ptr<EngngModel> InstanciateProblem(DataReader &dr, problemMode mode, int contextFlag, EngngModel *master = 0, bool parallelFlag = false);

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

    py::enum_<oofem::ValueModeType>(m, "ValueModeType")
        .value("VM_Unknown", oofem::ValueModeType::VM_Unknown)
        .value("VM_Total", oofem::ValueModeType::VM_Total)
        .value("VM_Velocity", oofem::ValueModeType::VM_Velocity)
        .value("VM_Acceleration", oofem::ValueModeType::VM_Acceleration)
        .value("VM_Incremental", oofem::ValueModeType::VM_Incremental)
        .value("VM_RhsTotal", oofem::ValueModeType::VM_RhsTotal)
        .value("VM_RhsIncremental", oofem::ValueModeType::VM_RhsIncremental)
        .value("VM_RhsInitial", oofem::ValueModeType::VM_RhsInitial)
        .value("VM_Intermediate", oofem::ValueModeType::VM_Intermediate)
        .value("VM_TotalIntrinsic", oofem::ValueModeType::VM_TotalIntrinsic)
    ;

    py::enum_<oofem::DofIDItem>(m, "DofIDItem")
      .value("Undef", oofem::DofIDItem::Undef)
      .value("D_u", oofem::DofIDItem::D_u)
      .value("D_v", oofem::DofIDItem::D_v)
      .value("D_w", oofem::DofIDItem::D_w)
      .value("R_u", oofem::DofIDItem::R_u)
      .value("R_v", oofem::DofIDItem::R_v)
      .value("R_w", oofem::DofIDItem::R_w)
      .value("V_u", oofem::DofIDItem::V_u)
      .value("V_v", oofem::DofIDItem::V_v)
      .value("V_w", oofem::DofIDItem::V_w)
      .value("T_f", oofem::DofIDItem::T_f)
      .value("P_f", oofem::DofIDItem::P_f)
      .value("G_0", oofem::DofIDItem::G_0)
      .value("G_1", oofem::DofIDItem::G_1)
      .value("C_1", oofem::DofIDItem::C_1)
      .value("W_u", oofem::DofIDItem::W_u)
      .value("W_v", oofem::DofIDItem::W_v)
      .value("W_w", oofem::DofIDItem::W_w)
      .value("Gamma", oofem::DofIDItem::Gamma)
      .value("D_u_edge_const", oofem::DofIDItem::D_u_edge_const)
      .value("D_u_edge_lin", oofem::DofIDItem::D_u_edge_lin)
      .value("D_v_edge_const", oofem::DofIDItem::D_v_edge_const)
      .value("D_v_edge_lin", oofem::DofIDItem::D_v_edge_lin)
      .value("Warp_PsiTheta", oofem::DofIDItem::Warp_PsiTheta)
      .value("Warp_Theta", oofem::DofIDItem::Warp_Theta)
      .value("LMP_u", oofem::DofIDItem::LMP_u)
      .value("LMP_v", oofem::DofIDItem::LMP_v)
      .value("LMP_w", oofem::DofIDItem::LMP_w)
      .value("Trac_u", oofem::DofIDItem::Trac_u)
      .value("Trac_v", oofem::DofIDItem::Trac_v)
      .value("Trac_w", oofem::DofIDItem::Trac_w)
      ;

    py::enum_<oofem::problemMode>(m, "problemMode")
        .value("processor", oofem::problemMode::_processor)
        .value("postProcessor", oofem::problemMode::_postProcessor)
    ;

}
