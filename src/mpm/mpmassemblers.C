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

#include "mpmassemblers.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "dictionary.h"
#include "verbose.h"
#include "classfactory.h"
#include "mathfem.h"
#include "assemblercallback.h"
#include "unknownnumberingscheme.h"
#include "dofdistributedprimaryfield.h"
#include "primaryfield.h"
#include "maskedprimaryfield.h"
#include "nrsolver.h"
#include "activebc.h"
#include "boundarycondition.h"
#include "boundaryload.h"
#include "outputmanager.h"
#include "mpm.h"

namespace oofem {

UPLhsAssembler :: UPLhsAssembler(double alpha, double deltaT) : 
    MatrixAssembler(), alpha(alpha), deltaT(deltaT)
{}


void UPLhsAssembler :: matrixFromElement(FloatMatrix &answer, Element &el, TimeStep *tStep) const
{
    FloatMatrix contrib;
    IntArray locu, locp;
    MPElement *e = dynamic_cast<MPElement*>(&el);
    int ndofs = e->giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();

    e->getLocalCodeNumbers (locu, Variable::VariableQuantity::Displacement);
    e->getLocalCodeNumbers (locp, Variable::VariableQuantity::Pressure);

    e->giveCharacteristicMatrix(contrib, MomentumBalance_StiffnessMatrix, tStep);
    contrib.times(this->alpha);
    answer.assemble(contrib, locu, locu);
    e->giveCharacteristicMatrix(contrib, MomentumBalance_PressureCouplingMatrix, tStep);
    contrib.times((-1.0)*this->alpha);
    answer.assemble(contrib, locu, locp);

    e->giveCharacteristicMatrix(contrib, MassBalance_PermeabilityMatrix, tStep);
    contrib.times((-1.0)*this->alpha*this->alpha*this->deltaT);
    answer.assemble(contrib, locp, locp);
    e->giveCharacteristicMatrix(contrib, MassBalance_CompresibilityMatrix, tStep);
    contrib.times((-1.0)*this->alpha);
    answer.assemble(contrib, locp, locp);
    e->giveCharacteristicMatrix(contrib, MassBalance_StressCouplingMatrix, tStep);
    contrib.times((-1.0)*this->alpha);
    answer.assemble(contrib, locp, locu);
}

void UPResidualAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    FloatArray contrib;
    IntArray locu, locp;
    MPElement *e = dynamic_cast<MPElement*>(&element);
    int ndofs = e->giveNumberOfDofs();
    vec.resize(ndofs);
    vec.zero();

    e->getLocalCodeNumbers (locu, Variable::VariableQuantity::Displacement);
    e->getLocalCodeNumbers (locp, Variable::VariableQuantity::Pressure);

    e->giveCharacteristicVector(contrib, MomentumBalance_StressResidual, mode, tStep);
    vec.assemble(contrib, locu);
    e->giveCharacteristicVector(contrib, MomentumBalance_PressureResidual, mode, tStep);
    contrib.times((-1.0));
    vec.assemble(contrib, locu);

    e->giveCharacteristicVector(contrib, MassBalance_StressRateResidual, mode, tStep);
    contrib.times(-1.0*alpha*deltaT);
    vec.assemble(contrib, locp);
    e->giveCharacteristicVector(contrib, MassBalance_PressureResidual, mode, tStep);
    contrib.times(-1.0*alpha*deltaT);
    vec.assemble(contrib, locp);
    e->giveCharacteristicVector(contrib, MassBalance_PressureRateResidual, mode, tStep);
    contrib.times(-1.0*alpha*deltaT);
    vec.assemble(contrib, locp);

    //vec.negated();
}


TMLhsAssembler :: TMLhsAssembler(double alpha, double deltaT) : 
    MatrixAssembler(), alpha(alpha), deltaT(deltaT)
{}


void TMLhsAssembler :: matrixFromElement(FloatMatrix &answer, Element &el, TimeStep *tStep) const
{
    FloatMatrix contrib;
    IntArray locu, loct;
    MPElement *e = dynamic_cast<MPElement*>(&el);
    int ndofs = e->giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();

    e->getLocalCodeNumbers (locu, Variable::VariableQuantity::Displacement);
    e->getLocalCodeNumbers (loct, Variable::VariableQuantity::Temperature);

    e->giveCharacteristicMatrix(contrib, MomentumBalance_StiffnessMatrix, tStep);
    contrib.times(this->alpha);
    answer.assemble(contrib, locu, locu);
    e->giveCharacteristicMatrix(contrib, MomentumBalance_ThermalCouplingMatrix, tStep);
    contrib.times(this->alpha);
    answer.assemble(contrib, locu, loct);

    e->giveCharacteristicMatrix(contrib, EnergyBalance_ConductivityMatrix, tStep);
    contrib.times(this->alpha);
    answer.assemble(contrib, loct, loct);
    e->giveCharacteristicMatrix(contrib, EnergyBalance_CapacityMatrix, tStep);
    contrib.times(1/tStep->giveTimeIncrement());
    answer.assemble(contrib, loct, loct);

    // @bp: experimental: evaluate element load linearization as part of residual
    // note: extend bctracker to track bc applied via sets and directly on elements!
    // node: ideally all external load assembled using term concept as part of residual vector! (including nodes)
    // this would allow for consistent assembly of linearized terms as well.
    BCTracker *bct = el.giveDomain()->giveBCTracker();
    BCTracker::entryListType bcList = bct->getElementRecords(el.giveNumber());
    // loop over all boundary conditions applied to the element
    for (BCTracker::entryListType::iterator it = bcList.begin(); it != bcList.end(); ++it) {
        BoundaryLoad *bc = dynamic_cast<BoundaryLoad*>(el.giveDomain()->giveBc((*it).bcNumber));
        int boundaryID = (*it).boundaryId;
        if (bc) {
            if (bc->giveType() == ConvectionBC) {
                e->giveCharacteristicMatrixFromBC(contrib, EnergyBalance_ConvectionBCMatrix, tStep, bc, boundaryID);
                if(contrib.isNotEmpty()) {
                    contrib.times(this->alpha);
                    if (bc->giveBCGeoType() == bcGeomType::SurfaceLoadBGT) {
                        e->getSurfaceElementCodeNumbers(loct, Variable::VariableQuantity::Temperature, boundaryID);
                    } else {
                        e->getEdgeElementCodeNumbers(loct, Variable::VariableQuantity::Temperature, boundaryID);
                    }
                    answer.assemble(contrib, loct, loct);
                }
            }
        }
    
    }

}

void TMResidualAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    FloatArray contrib;
    IntArray locu, loct;
    MPElement *e = dynamic_cast<MPElement*>(&element);
    int ndofs = e->giveNumberOfDofs();
    vec.resize(ndofs);
    vec.zero();

    e->getLocalCodeNumbers (locu, Variable::VariableQuantity::Displacement);
    e->getLocalCodeNumbers (loct, Variable::VariableQuantity::Temperature);

    e->giveCharacteristicVector(contrib, MomentumBalance_StressResidual, mode, tStep);
    vec.assemble(contrib, locu);
    
    e->giveCharacteristicVector(contrib, EnergyBalance_Residual, mode, tStep);
    vec.assemble(contrib, loct);

    BCTracker *bct = e->giveDomain()->giveBCTracker();
    BCTracker::entryListType bcList = bct->getElementRecords(element.giveNumber());
    // loop over all boundary conditions applied to the element
    for (BCTracker::entryListType::iterator it = bcList.begin(); it != bcList.end(); ++it) {
        BoundaryLoad *bc = dynamic_cast<BoundaryLoad*>(element.giveDomain()->giveBc((*it).bcNumber));
        int boundaryID = (*it).boundaryId;
        if (bc) {
            if (bc->giveType() == ConvectionBC) {
                e->giveCharacteristicVectorFromBC(contrib, EnergyBalance_ConvectionBCResidual, mode, tStep, bc, boundaryID);
                if(contrib.isNotEmpty()) {
                    if (bc->giveBCGeoType() == bcGeomType::SurfaceLoadBGT) {
                        e->getSurfaceElementCodeNumbers(loct, Variable::VariableQuantity::Temperature, boundaryID);
                    } else {
                        e->getEdgeElementCodeNumbers(loct, Variable::VariableQuantity::Temperature, boundaryID);
                    }
                    vec.assemble(contrib, loct);
                }
            }
        }
    
    }

}

} // end namespace oofem
