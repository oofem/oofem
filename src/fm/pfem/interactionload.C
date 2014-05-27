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

#include "interactionload.h"
#include "floatarray.h"
#include "mathfem.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "timestep.h"
#include "fluidstructureproblem.h"
#include "pfem.h"
#include "loadtimefunction.h"

namespace oofem {

REGISTER_BoundaryCondition( InteractionLoad );

IRResultType
InteractionLoad :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    LinearEdgeLoad :: initializeFrom(ir);

	IR_GIVE_FIELD(ir, coupledParticles, _IFT_InteractionLoad_CoupledParticles);
 ///   if ( componentArray.giveSize() != nDofs * 2 ) {
 ///       _error("instanciateFrom: componentArray size mismatch");
 ///   }
 ///
 ///   int fType = 0;
 ///   IR_GIVE_OPTIONAL_FIELD(ir, fType, _IFT_LinearEdgeLoad_formulation);
 ///   if ( fType == 1 ) {
 ///       this->formulation = FT_Global;
 ///       // read start and end coordinates
 ///       IR_GIVE_FIELD(ir, startCoords, _IFT_LinearEdgeLoad_startcoord);
 ///       IR_GIVE_FIELD(ir, endCoords, _IFT_LinearEdgeLoad_endcoord);
 ///       if ( startCoords.isEmpty() || endCoords.isEmpty() ) {
 ///           _error("instanciateFrom: coordinates not specified");
 ///       }
 ///   } else {
 ///       this->formulation = FT_Entity;
 ///   }

    return IRRT_OK;
}


void InteractionLoad :: giveInputRecord(DynamicInputRecord& input)
{
    LinearEdgeLoad :: giveInputRecord ( input );
 ///   input.setField(this->formulation, _IFT_LinearEdgeLoad_formulation);
 ///   if ( this->formulation == FT_Global ) {
 ///       input.setField(this->startCoords, _IFT_LinearEdgeLoad_startcoord);
 ///       input.setField(this->endCoords, _IFT_LinearEdgeLoad_endcoord);
 ///   }
}

void
InteractionLoad::computeValueAt(FloatArray &answer, TimeStep *tStep, FloatArray &coords, ValueModeType mode)
{
	FloatArray pressureArray(this->componentArray);
	//tStep->giveEngngModel();
	FluidStructureProblem *fsiProblem = dynamic_cast<FluidStructureProblem*>(domain->giveEngngModel()->giveMasterEngngModel());
	if (fsiProblem) {
		for ( int i = 1; i <= fsiProblem->giveNumberOfSlaveProblems(); i++ ) {
			PFEM *pfem = dynamic_cast<PFEM*>(fsiProblem->giveSlaveProblem(i));
			if (pfem) {
				for ( int j = 1; j <= coupledParticles.giveSize(); j++) {
					DofManager *dman = pfem->giveDomain(1)->giveDofManager(coupledParticles.at(j));
					Dof *pressureDof = dman->giveDofWithID(P_f);
					double pressureValue = pfem->giveUnknownComponent(VM_Total, tStep, pfem->giveDomain(1), pressureDof);
					pressureValue = pressureValue > 0 ? pressureValue : 0.0;
					for ( int k = 1; k <= nDofs; k++) {
						pressureArray.at(nDofs * (j - 1) + k) *= pressureValue;
					}
				}
			}
		}
	}
	// Evaluates the value at specific integration point
    int i, j, nSize;
    double value, factor;
    FloatArray N;

    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        _error("computeValueAt: unknown mode");
    }

    answer.resize(this->nDofs);

    this->computeNArray(N, coords);
    nSize = N.giveSize();

    if ( ( pressureArray.giveSize() / nSize ) != nDofs ) {
        _error("computeValueAt: componentArray size mismatch");
    }

    for ( i = 1; i <= nDofs; i++ ) {
        for ( value = 0., j = 1; j <= nSize; j++ ) {
            value += N.at(j) * pressureArray.at(i + ( j - 1 ) * nDofs);
        }

        answer.at(i) = value;
    }

   
    factor = this->giveLoadTimeFunction()->evaluate(tStep, mode);

    answer.times(factor);
}

void
InteractionLoad :: computeNArray(FloatArray &answer, FloatArray &coords) const
{
	LinearEdgeLoad :: computeNArray(answer, coords);
 ///   // compute local isoparametric coordinates of given point
 ///   double ksi;
 ///
 ///   if ( formulation == FT_Global ) {
 ///       int i;
 ///       double length = endCoords.distance(startCoords);
 ///       double dl     = coords.distance(startCoords);
 ///       double eta = dl / length;
 ///       ksi    = ( dl - 0.5 * length ) / ( 0.5 * length );
 ///       FloatArray dir = endCoords;
 ///
 ///       dir.subtract(startCoords);
 ///
 ///       if ( ( ksi < -1.0 ) ||  ( ksi > 1.0 ) ) {
 ///           _warning2("computeNArray: point out of receiver, skipped", 1);
 ///           answer.resize(2);
 ///           answer.zero();
 ///       }
 ///
 ///       for ( i = 1; i <= dir.giveSize(); i++ ) {
 ///           if ( fabs( startCoords.at(i) + dir.at(i) * eta - coords.at(i) ) > 1.e-6 ) {
 ///               _warning2("computeNArray: point out of receiver, skipped", 1);
 ///               answer.resize(2);
 ///               answer.zero();
 ///           }
 ///       }
 ///   } else {
 ///       ksi = coords.at(1);
 ///   }
 ///
 ///   double n1, n2;
 ///
 ///   n1  = ( 1. - ksi ) * 0.5;
 ///   n2  = ( 1. + ksi ) * 0.5;
 ///
 ///   answer.resize(2);
 ///
 ///   answer.at(1) = n1;
 ///   answer.at(2) = n2;
}
} // end namespace oofem
