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
#include "function.h"

namespace oofem {
REGISTER_BoundaryCondition(InteractionLoad);

IRResultType
InteractionLoad :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, coupledParticles, _IFT_InteractionLoad_CoupledParticles);


    return LinearEdgeLoad :: initializeFrom(ir);
}


void InteractionLoad :: giveInputRecord(DynamicInputRecord &input)
{
    LinearEdgeLoad :: giveInputRecord(input);
}

void
InteractionLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
{
    int nSize, nDofs;
    double factor;
    FloatArray N;

    // Evaluates the value at specific integration point
    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        OOFEM_ERROR("computeValueAt: unknown mode");
    }

    this->computeNArray(N, coords);
    nSize = N.giveSize();

    nDofs = this->componentArray.giveSize() / nSize;

    FloatArray pressureArray = this->componentArray;
    FluidStructureProblem *fsiProblem = dynamic_cast< FluidStructureProblem * >( domain->giveEngngModel()->giveMasterEngngModel() );
    if ( fsiProblem ) {
        for ( int i = 1; i <= fsiProblem->giveNumberOfSlaveProblems(); i++ ) {
            PFEM *pfem = dynamic_cast< PFEM * >( fsiProblem->giveSlaveProblem(i) );
            if ( pfem ) {
                for ( int j = 1; j <= coupledParticles.giveSize(); j++ ) {
                    DofManager *dman = pfem->giveDomain(1)->giveDofManager( coupledParticles.at(j) );
                    Dof *pressureDof = dman->giveDofWithID(P_f);
                    double pressureValue = pfem->giveUnknownComponent(VM_Total, tStep, pfem->giveDomain(1), pressureDof);
                    pressureValue = pressureValue > 0 ? pressureValue : 0.0;
                    for ( int k = 1; k <= nDofs; k++ ) {
                        pressureArray.at(nDofs * ( j - 1 ) + k) *= pressureValue;
                    }
                }
            }
        }
    }

    answer.resize(nDofs);

    if ( ( pressureArray.giveSize() / nSize ) != nDofs ) {
        OOFEM_ERROR("computeValueAt: componentArray size mismatch");
    }

    for ( int i = 1; i <= nDofs; i++ ) {
        double value = 0.;
        for ( int j = 1; j <= nSize; j++ ) {
            value += N.at(j) * pressureArray.at(i + ( j - 1 ) * nDofs);
        }

        answer.at(i) = value;
    }

    factor = this->giveTimeFunction()->evaluate(tStep, mode);

    answer.times(factor);
}

void
InteractionLoad :: computeNArray(FloatArray &answer, const FloatArray &coords) const
{
    LinearEdgeLoad :: computeNArray(answer, coords);
}
} // end namespace oofem
