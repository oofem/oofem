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
#include "chartype.h"
#include "elementgeometrytype.h"
#include "internalstatetype.h"

#include "integrationrule.h"
#include "gausspoint.h"
#include "inputrecord.h"
#include "dynamicinputrecord.h"

#include "classfactory.h"
#include "unknownnumberingscheme.h"


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

    py::class_<oofem::InputRecord>(m, "InputRecord")
    ;

    typedef const char *InputFieldType;
    py::class_<oofem::DynamicInputRecord, oofem::InputRecord>(m, "DynamicInputRecord")
        .def(py::init<std::string, int>(), py::arg("answer") = "", py::arg("value")=0)
        .def("finish", &oofem::DynamicInputRecord::finish, py::arg("wrn")=true)
        .def("setRecordKeywordField", &oofem::DynamicInputRecord::setRecordKeywordField)
        .def("setRecordKeywordNumber", &oofem::DynamicInputRecord::setRecordKeywordNumber)
        .def("setField", (void (oofem::DynamicInputRecord::*)(int, InputFieldType)) &oofem::DynamicInputRecord::setField) 
        .def("setField", (void (oofem::DynamicInputRecord::*)(double, InputFieldType)) &oofem::DynamicInputRecord::setField) 
        .def("setField", (void (oofem::DynamicInputRecord::*)(bool, InputFieldType)) &oofem::DynamicInputRecord::setField) 
        .def("setField", (void (oofem::DynamicInputRecord::*)(std::string, InputFieldType)) &oofem::DynamicInputRecord::setField) 
        .def("setField", (void (oofem::DynamicInputRecord::*)(oofem::FloatArray, InputFieldType)) &oofem::DynamicInputRecord::setField) 
        .def("setField", (void (oofem::DynamicInputRecord::*)(oofem::FloatMatrix, InputFieldType)) &oofem::DynamicInputRecord::setField) 
        .def("setField", (void (oofem::DynamicInputRecord::*)(oofem::FloatMatrix, InputFieldType)) &oofem::DynamicInputRecord::setField) 
        .def("setField", (void (oofem::DynamicInputRecord::*)(oofem::IntArray, InputFieldType)) &oofem::DynamicInputRecord::setField) 
        .def("setField", (void (oofem::DynamicInputRecord::*)(std::vector<std::string>, InputFieldType)) &oofem::DynamicInputRecord::setField) 
        .def("setField", (void (oofem::DynamicInputRecord::*)(InputFieldType)) &oofem::DynamicInputRecord::setField) 
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
        .def("hasCoordinates", &oofem::DofManager::hasCoordinates)
        .def("giveCoordinates", &oofem::DofManager::giveCoordinates, py::return_value_policy::reference)
        .def("appendDof", &oofem::DofManager::appendDof, py::keep_alive<1, 2>())
        .def("removeDof", &oofem::DofManager::removeDof)
        .def("hasDofID", &oofem::DofManager::hasDofID)
    ;

    py::class_<oofem::Element, oofem::FEMComponent>(m, "Element")
        .def("giveLocationArray", (void (oofem::Element::*)(oofem::IntArray &, const oofem::UnknownNumberingScheme &, oofem::IntArray *dofIds) const) &oofem::Element::giveLocationArray)
        .def("giveLocationArray", (void (oofem::Element::*)(oofem::IntArray &, const oofem::IntArray &, const oofem::UnknownNumberingScheme &, oofem::IntArray *) const) &oofem::Element::giveLocationArray)
        .def("giveCharacteristicMatrix", &oofem::Element::giveCharacteristicMatrix)
        .def("giveCharacteristicVector", &oofem::Element::giveCharacteristicVector)
        .def("giveCharacteristicValue", &oofem::Element::giveCharacteristicValue)
        .def("computeVectorOf", (void (oofem::Element::*)(oofem::ValueModeType , oofem::TimeStep *, oofem::FloatArray &)) &oofem::Element::computeVectorOf)
        .def("computeVectorOf", (void (oofem::Element::*)(const oofem::IntArray &, oofem::ValueModeType , oofem::TimeStep *, oofem::FloatArray &, bool)) &oofem::Element::computeVectorOf)
        .def("computeVectorOfPrescribed", (void (oofem::Element::*)(oofem::ValueModeType , oofem::TimeStep *, oofem::FloatArray &)) &oofem::Element::computeVectorOfPrescribed)
        .def("computeVectorOfPrescribed", (void (oofem::Element::*)(const oofem::IntArray &, oofem::ValueModeType , oofem::TimeStep *, oofem::FloatArray & )) &oofem::Element::computeVectorOfPrescribed)
        .def("giveDofManagerNumber", &oofem::Element::giveDofManagerNumber)
        .def("giveDofManArray", &oofem::Element::giveDofManArray, py::return_value_policy::reference)
        .def("giveDofManager", &oofem::Element::giveDofManager, py::return_value_policy::reference)
        .def("giveMaterial", &oofem::Element::giveMaterial, py::return_value_policy::reference)
        .def("giveCrossSection", &oofem::Element::giveCrossSection, py::return_value_policy::reference)
        .def("giveMaterialNumber", &oofem::Element::giveMaterialNumber)
        .def("giveNumberOfDofManagers", &oofem::Element::giveNumberOfDofManagers)
        .def("giveNumberOfNodes", &oofem::Element::giveNumberOfNodes)
        .def("giveRegionNumber", &oofem::Element::giveRegionNumber)
        .def("updateYourself", &oofem::Element::updateYourself)
        .def("isActivated", &oofem::Element::isActivated)
        .def("isCast", &oofem::Element::isCast)
        .def("giveGeometryType", &oofem::Element::giveGeometryType)
        .def("giveIntegrationRule", &oofem::Element::giveIntegrationRule, py::return_value_policy::reference)
        .def("giveDefaultIntegrationRulePtr", &oofem::Element::giveDefaultIntegrationRulePtr, py::return_value_policy::reference)
        .def("giveIPValue", &oofem::Element::giveIPValue)
        .def("giveLabel", &oofem::Element::giveLabel)
        .def("initializeFrom", &oofem::Element::initializeFrom)
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

    py::class_<oofem::UnknownNumberingScheme>(m, "UnknownNumberingScheme")
    ;

    py::class_<oofem::IntegrationRule>(m, "IntegrationRule")
    ;
    py::class_<oofem::GaussPoint>(m, "GaussPoint")
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

    py::enum_<oofem::CharType>(m, "CharType")
      .value("UnknownCharType", oofem::CharType::UnknownCharType)
      .value("StiffnessMatrix", oofem::CharType::StiffnessMatrix)
      .value("TangentStiffnessMatrix", oofem::CharType::TangentStiffnessMatrix)
      .value("SecantStiffnessMatrix", oofem::CharType::SecantStiffnessMatrix)
      .value("ElasticStiffnessMatrix", oofem::CharType::ElasticStiffnessMatrix)
      .value("MassMatrix", oofem::CharType::MassMatrix)
      .value("LumpedMassMatrix", oofem::CharType::LumpedMassMatrix)
      .value("ConductivityMatrix", oofem::CharType::ConductivityMatrix)
      .value("CapacityMatrix", oofem::CharType::CapacityMatrix)
      .value("InitialStressMatrix", oofem::CharType::InitialStressMatrix)
      .value("ExternalForcesVector", oofem::CharType::ExternalForcesVector)
      .value("InternalForcesVector", oofem::CharType::InternalForcesVector)
      .value("LastEquilibratedInternalForcesVector", oofem::CharType::LastEquilibratedInternalForcesVector)
      .value("InertiaForcesVector", oofem::CharType::InertiaForcesVector)
      .value("AuxVelocityLhs", oofem::CharType::AuxVelocityLhs)
      .value("VelocityLhs", oofem::CharType::VelocityLhs)
      .value("PressureGradientMatrix", oofem::CharType::PressureGradientMatrix)
      .value("DivergenceMatrix", oofem::CharType::DivergenceMatrix)
      .value("VelocityLaplacianMatrix", oofem::CharType::VelocityLaplacianMatrix)
      .value("PressureLaplacianMatrix", oofem::CharType::PressureLaplacianMatrix)
      .value("StabilizationMassMatrix", oofem::CharType::StabilizationMassMatrix)
      .value("PressureGradientVector", oofem::CharType::PressureGradientVector)
      .value("MassVelocityVector", oofem::CharType::MassVelocityVector)
      .value("MassAuxVelocityVector", oofem::CharType::MassAuxVelocityVector)
      .value("LaplacePressureVector", oofem::CharType::LaplacePressureVector)
      .value("LaplaceVelocityVector", oofem::CharType::LaplaceVelocityVector)
      .value("DivergenceAuxVelocityVector", oofem::CharType::DivergenceAuxVelocityVector)
      .value("DivergenceVelocityVector", oofem::CharType::DivergenceVelocityVector)
      ;

    py::enum_<oofem::Element_Geometry_Type>(m, "Element_Geometry_Type")
      .value("EGT_point", oofem::Element_Geometry_Type::EGT_point)
      .value("EGT_line_1", oofem::Element_Geometry_Type::EGT_line_1)
      .value("EGT_line_2", oofem::Element_Geometry_Type::EGT_line_2)
      .value("EGT_triangle_1", oofem::Element_Geometry_Type::EGT_triangle_1)
      .value("EGT_triangle_2", oofem::Element_Geometry_Type::EGT_triangle_2)
      .value("EGT_quad_1", oofem::Element_Geometry_Type::EGT_quad_1)
      .value("EGT_quad_1_interface", oofem::Element_Geometry_Type::EGT_quad_1_interface)
      .value("EGT_quad_21_interface", oofem::Element_Geometry_Type::EGT_quad_21_interface)
      .value("EGT_quad_2", oofem::Element_Geometry_Type::EGT_quad_2)
      .value("EGT_quad9_2", oofem::Element_Geometry_Type::EGT_quad9_2)
      .value("EGT_tetra_1", oofem::Element_Geometry_Type::EGT_tetra_1)
      .value("EGT_tetra_2", oofem::Element_Geometry_Type::EGT_tetra_2)
      .value("EGT_hexa_1", oofem::Element_Geometry_Type::EGT_hexa_1)
      .value("EGT_hexa_2", oofem::Element_Geometry_Type::EGT_hexa_2)
      .value("EGT_hexa_27", oofem::Element_Geometry_Type::EGT_hexa_27)
      .value("EGT_wedge_1", oofem::Element_Geometry_Type::EGT_wedge_1)
      .value("EGT_wedge_2", oofem::Element_Geometry_Type::EGT_wedge_2)
      .value("EGT_Composite", oofem::Element_Geometry_Type::EGT_Composite)
      .value("EGT_unknown", oofem::Element_Geometry_Type::EGT_unknown)
      ;

    py::enum_<oofem::InternalStateType>(m, "InternalStateType")
      .value("IST_Undefined", oofem::InternalStateType::IST_Undefined)
      .value("IST_StressTensor", oofem::InternalStateType::IST_StressTensor)
      .value("IST_PrincipalStressTensor", oofem::InternalStateType::IST_PrincipalStressTensor)
      .value("IST_PrincipalStressTempTensor", oofem::InternalStateType::IST_PrincipalStressTempTensor)
      .value("IST_StrainTensor", oofem::InternalStateType::IST_StrainTensor)
      .value("IST_PrincipalStrainTensor", oofem::InternalStateType::IST_PrincipalStrainTensor)
      .value("IST_PrincipalStrainTempTensor", oofem::InternalStateType::IST_PrincipalStrainTempTensor)
      .value("IST_BeamForceMomentTensor", oofem::InternalStateType::IST_BeamForceMomentTensor)
      .value("IST_BeamStrainCurvatureTensor", oofem::InternalStateType::IST_BeamStrainCurvatureTensor)
      .value("IST_ShellMomentTensor", oofem::InternalStateType::IST_ShellMomentTensor)
      .value("IST_ShellForceTensor", oofem::InternalStateType::IST_ShellForceTensor)
      .value("IST_CurvatureTensor", oofem::InternalStateType::IST_CurvatureTensor)
      .value("IST_DisplacementVector", oofem::InternalStateType::IST_DisplacementVector)
      .value("IST_DamageTensor", oofem::InternalStateType::IST_DamageTensor)
      .value("IST_DamageInvTensor", oofem::InternalStateType::IST_DamageInvTensor)
      .value("IST_PrincipalDamageTensor", oofem::InternalStateType::IST_PrincipalDamageTensor)
      .value("IST_PrincipalDamageTempTensor", oofem::InternalStateType::IST_PrincipalDamageTempTensor)
      .value("IST_CrackState", oofem::InternalStateType::IST_CrackState)
      .value("IST_StressTensorTemp", oofem::InternalStateType::IST_StressTensorTemp)
      .value("IST_StrainTensorTemp", oofem::InternalStateType::IST_StrainTensorTemp)
      .value("IST_ShellForceTensorTemp", oofem::InternalStateType::IST_ShellForceTensorTemp)
      .value("IST_ShellMomentTensorTemp", oofem::InternalStateType::IST_ShellMomentTensorTemp)
      .value("IST_CurvatureTensorTemp", oofem::InternalStateType::IST_CurvatureTensorTemp)
      .value("IST_DisplacementVectorTemp", oofem::InternalStateType::IST_DisplacementVectorTemp)
      .value("IST_DamageTensorTemp", oofem::InternalStateType::IST_DamageTensorTemp)
      .value("IST_DamageInvTensorTemp", oofem::InternalStateType::IST_DamageInvTensorTemp)
      .value("IST_CrackStateTemp", oofem::InternalStateType::IST_CrackStateTemp)
      .value("IST_PlasticStrainTensor", oofem::InternalStateType::IST_PlasticStrainTensor)
      .value("IST_PrincipalPlasticStrainTensor", oofem::InternalStateType::IST_PrincipalPlasticStrainTensor)
      .value("IST_CylindricalStressTensor", oofem::InternalStateType::IST_CylindricalStressTensor)
      .value("IST_CylindricalStrainTensor", oofem::InternalStateType::IST_CylindricalStrainTensor)
      .value("IST_MaxEquivalentStrainLevel", oofem::InternalStateType::IST_MaxEquivalentStrainLevel)
      .value("IST_ErrorIndicatorLevel", oofem::InternalStateType::IST_ErrorIndicatorLevel)
      .value("IST_InternalStressError", oofem::InternalStateType::IST_InternalStressError)
      .value("IST_PrimaryUnknownError", oofem::InternalStateType::IST_PrimaryUnknownError)
      .value("IST_RelMeshDensity", oofem::InternalStateType::IST_RelMeshDensity)
      .value("IST_MicroplaneDamageValues", oofem::InternalStateType::IST_MicroplaneDamageValues)
      .value("IST_Temperature", oofem::InternalStateType::IST_Temperature)
      .value("IST_MassConcentration_1", oofem::InternalStateType::IST_MassConcentration_1)
      .value("IST_HydrationDegree", oofem::InternalStateType::IST_HydrationDegree)
      .value("IST_Humidity", oofem::InternalStateType::IST_Humidity)
      .value("IST_Velocity", oofem::InternalStateType::IST_Velocity)
      .value("IST_Pressure", oofem::InternalStateType::IST_Pressure)
      .value("IST_VOFFraction", oofem::InternalStateType::IST_VOFFraction)
      .value("IST_Density", oofem::InternalStateType::IST_Density)
      .value("IST_MaterialInterfaceVal", oofem::InternalStateType::IST_MaterialInterfaceVal)
      .value("IST_MaterialNumber", oofem::InternalStateType::IST_MaterialNumber)
      .value("IST_ElementNumber", oofem::InternalStateType::IST_ElementNumber)
      .value("IST_BoneVolumeFraction", oofem::InternalStateType::IST_BoneVolumeFraction)
      .value("IST_PlasStrainEnerDens", oofem::InternalStateType::IST_PlasStrainEnerDens)
      .value("IST_ElasStrainEnerDens", oofem::InternalStateType::IST_ElasStrainEnerDens)
      .value("IST_TotalStrainEnerDens", oofem::InternalStateType::IST_TotalStrainEnerDens)
      .value("IST_DamageScalar", oofem::InternalStateType::IST_DamageScalar)
      .value("IST_MaterialOrientation_x", oofem::InternalStateType::IST_MaterialOrientation_x)
      .value("IST_MaterialOrientation_y", oofem::InternalStateType::IST_MaterialOrientation_y)
      .value("IST_MaterialOrientation_z", oofem::InternalStateType::IST_MaterialOrientation_z)
      .value("IST_TemperatureFlow", oofem::InternalStateType::IST_TemperatureFlow)
      .value("IST_MassConcentrationFlow_1", oofem::InternalStateType::IST_MassConcentrationFlow_1)
      .value("IST_HumidityFlow", oofem::InternalStateType::IST_HumidityFlow)
      .value("IST_CrackStatuses", oofem::InternalStateType::IST_CrackStatuses)
      .value("IST_CrackedFlag", oofem::InternalStateType::IST_CrackedFlag)
      .value("IST_CrackDirs", oofem::InternalStateType::IST_CrackDirs)
      .value("IST_CumPlasticStrain", oofem::InternalStateType::IST_CumPlasticStrain)
      .value("IST_CumPlasticStrain_2", oofem::InternalStateType::IST_CumPlasticStrain_2)
      .value("IST_StressWorkDensity", oofem::InternalStateType::IST_StressWorkDensity)
      .value("IST_DissWorkDensity", oofem::InternalStateType::IST_DissWorkDensity)
      .value("IST_FreeEnergyDensity", oofem::InternalStateType::IST_FreeEnergyDensity)
      .value("IST_ThermalConductivityIsotropic", oofem::InternalStateType::IST_ThermalConductivityIsotropic)
      .value("IST_HeatCapacity", oofem::InternalStateType::IST_HeatCapacity)
      .value("IST_AverageTemperature", oofem::InternalStateType::IST_AverageTemperature)
      .value("IST_YoungModulusVirginPaste", oofem::InternalStateType::IST_YoungModulusVirginPaste)
      .value("IST_PoissonRatioVirginPaste", oofem::InternalStateType::IST_PoissonRatioVirginPaste)
      .value("IST_YoungModulusConcrete", oofem::InternalStateType::IST_YoungModulusConcrete)
      .value("IST_PoissonRatioConcrete", oofem::InternalStateType::IST_PoissonRatioConcrete)
      .value("IST_VolumetricPlasticStrain", oofem::InternalStateType::IST_VolumetricPlasticStrain)
      .value("IST_DeviatoricStrain", oofem::InternalStateType::IST_DeviatoricStrain)
      .value("IST_DeviatoricStress", oofem::InternalStateType::IST_DeviatoricStress)
      .value("IST_Viscosity", oofem::InternalStateType::IST_Viscosity)
      .value("IST_CharacteristicLength", oofem::InternalStateType::IST_CharacteristicLength)
      .value("IST_DeviatoricStrainMeasure", oofem::InternalStateType::IST_DeviatoricStrainMeasure)
      .value("IST_DeviatoricStressMeasure", oofem::InternalStateType::IST_DeviatoricStressMeasure)
      .value("IST_vonMisesStress", oofem::InternalStateType::IST_vonMisesStress)
      .value("IST_CrackVector", oofem::InternalStateType::IST_CrackVector)
      .value("IST_PressureGradient", oofem::InternalStateType::IST_PressureGradient)
      .value("IST_DissWork", oofem::InternalStateType::IST_DissWork)
      .value("IST_DeltaDissWork", oofem::InternalStateType::IST_DeltaDissWork)
      .value("IST_StressCapPos", oofem::InternalStateType::IST_StressCapPos)
      .value("IST_TangentNorm", oofem::InternalStateType::IST_TangentNorm)
      .value("IST_Tangent", oofem::InternalStateType::IST_Tangent)
      .value("IST_DirectorField", oofem::InternalStateType::IST_DirectorField)
      .value("IST_CrackWidth", oofem::InternalStateType::IST_CrackWidth)
      .value("IST_DeformationGradientTensor", oofem::InternalStateType::IST_DeformationGradientTensor)
      .value("IST_FirstPKStressTensor", oofem::InternalStateType::IST_FirstPKStressTensor)
      .value("IST_XFEMEnrichment", oofem::InternalStateType::IST_XFEMEnrichment)
      .value("IST_XFEMNumIntersecPoints", oofem::InternalStateType::IST_XFEMNumIntersecPoints)
      .value("IST_XFEMLevelSetPhi", oofem::InternalStateType::IST_XFEMLevelSetPhi)
      .value("IST_Maturity", oofem::InternalStateType::IST_Maturity)
      .value("IST_CauchyStressTensor", oofem::InternalStateType::IST_CauchyStressTensor)
      .value("IST_InterfaceJump", oofem::InternalStateType::IST_InterfaceJump)
      .value("IST_InterfaceTraction", oofem::InternalStateType::IST_InterfaceTraction)
      .value("IST_InterfaceFirstPKTraction", oofem::InternalStateType::IST_InterfaceFirstPKTraction)
      .value("IST_StressTensor_Reduced", oofem::InternalStateType::IST_StressTensor_Reduced)
      .value("IST_StrainTensor_Reduced", oofem::InternalStateType::IST_StrainTensor_Reduced)
      .value("IST_CrossSectionNumber", oofem::InternalStateType::IST_CrossSectionNumber)
      .value("IST_ShellStrainTensor", oofem::InternalStateType::IST_ShellStrainTensor)
      .value("IST_AbaqusStateVector", oofem::InternalStateType::IST_AbaqusStateVector)
      .value("IST_AutogenousShrinkageTensor", oofem::InternalStateType::IST_AutogenousShrinkageTensor)
      .value("IST_DryingShrinkageTensor", oofem::InternalStateType::IST_DryingShrinkageTensor)
      .value("IST_TotalShrinkageTensor", oofem::InternalStateType::IST_TotalShrinkageTensor)
      .value("IST_ThermalStrainTensor", oofem::InternalStateType::IST_ThermalStrainTensor)
      .value("IST_CreepStrainTensor", oofem::InternalStateType::IST_CreepStrainTensor)
      .value("IST_TensileStrength", oofem::InternalStateType::IST_TensileStrength)
      .value("IST_ResidualTensileStrength", oofem::InternalStateType::IST_ResidualTensileStrength)
      .value("IST_LocalEquivalentStrain", oofem::InternalStateType::IST_LocalEquivalentStrain)
      .value("IST_CrackIndex", oofem::InternalStateType::IST_CrackIndex)
      .value("IST_EigenStrainTensor", oofem::InternalStateType::IST_EigenStrainTensor)
      .value("IST_CrackStrainTensor", oofem::InternalStateType::IST_CrackStrainTensor)
      .value("IST_2ndCrackWidth", oofem::InternalStateType::IST_2ndCrackWidth)
      .value("IST_2ndCrackVector", oofem::InternalStateType::IST_2ndCrackVector)
      .value("IST_3rdCrackWidth", oofem::InternalStateType::IST_3rdCrackWidth)
      .value("IST_3rdCrackVector", oofem::InternalStateType::IST_3rdCrackVector)
      .value("IST_FiberStressLocal", oofem::InternalStateType::IST_FiberStressLocal)
      .value("IST_FiberStressNL", oofem::InternalStateType::IST_FiberStressNL)
      .value("IST_EnergyMassCapacity", oofem::InternalStateType::IST_EnergyMassCapacity)
      .value("IST_PrincStressVector1", oofem::InternalStateType::IST_PrincStressVector1)
      .value("IST_PrincStressVector2", oofem::InternalStateType::IST_PrincStressVector2)
      .value("IST_PrincStressVector3", oofem::InternalStateType::IST_PrincStressVector3)
      .value("IST_InterfaceNormal", oofem::InternalStateType::IST_InterfaceNormal)
      .value("IST_MomentTensor", oofem::InternalStateType::IST_MomentTensor)
      .value("IST_MomentTensorTemp", oofem::InternalStateType::IST_MomentTensorTemp)
      .value("IST_YieldStrength", oofem::InternalStateType::IST_YieldStrength)
      .value("IST_ElasticStrainTensor", oofem::InternalStateType::IST_ElasticStrainTensor)
      .value("IST_MoistureContent", oofem::InternalStateType::IST_MoistureContent)
      .value("IST_CrackStatusesTemp", oofem::InternalStateType::IST_CrackStatusesTemp)
      .value("IST_CrackSlip", oofem::InternalStateType::IST_CrackSlip)
      .value("IST_EquivalentTime", oofem::InternalStateType::IST_EquivalentTime)
      .value("IST_IncrementCreepModulus", oofem::InternalStateType::IST_IncrementCreepModulus)
      ;
}
