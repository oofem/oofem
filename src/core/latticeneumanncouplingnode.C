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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "latticeneumanncouplingnode.h"
#include "dof.h"
#include "slavedof.h"
#include "simpleslavedof.h"
#include "nodalload.h"
#include "timestep.h"

#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "verbose.h"
#include "datastream.h"
#include "contextioerr.h"
#include "staggeredproblem.h"

#ifdef __SM_MODULE
 #include "../sm/Elements/LatticeElements/latticestructuralelement.h"
#endif

#include <math.h>
#include <stdlib.h>
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_DofManager(LatticeNeumannCouplingNode);

LatticeNeumannCouplingNode :: LatticeNeumannCouplingNode(int n, Domain *aDomain) :
    Node(n, aDomain), directionVector(), couplingNodes()
{}

LatticeNeumannCouplingNode :: ~LatticeNeumannCouplingNode()
// Destructor.
{}

void
LatticeNeumannCouplingNode :: initializeFrom(InputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
    FloatArray triplets;

    Node :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, directionVector, _IFT_LatticeNeumannCouplingNode_direction);
    IR_GIVE_FIELD(ir, couplingNodes, _IFT_LatticeNeumannCouplingNode_couplingnodes);
}


void
LatticeNeumannCouplingNode :: computeLoadVectorAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Computes the vector of the nodal loads of the receiver.
{
    FloatArray contribution;

    if ( this->giveLoadArray()->isEmpty() ) {
        answer.resize(0);
    } else {
        answer.resize(0);
        int nLoads = loadArray.giveSize();     // the node may be subjected
        for ( int i = 1; i <= nLoads; i++ ) {   // to more than one load
            int n = loadArray.at(i);
            Load *loadN = domain->giveLoad(n);
            computeLoadVector(contribution, loadN, ExternalForcesVector, stepN, mode);
            answer.add(contribution);
        }
    }

    computeLoadCouplingContribution(contribution, stepN);
    answer.add(contribution);
}

void LatticeNeumannCouplingNode :: printYourself()
// Prints the receiver on screen.
{
    int i;
    double x, y;

    x = this->giveCoordinate(1);
    y = this->giveCoordinate(2);
    printf("LatticeNeumannCouplingNode %d    coord : x %f  y %f\n", number, x, y);
    for ( i = 0; i < this->giveNumberOfDofs(); i++ ) {
        if ( dofArray [ i ] ) {
            dofArray [ i ]->printYourself();
        } else {
            printf("dof %d is nil \n", i + 1);
        }
    }

    loadArray.printYourself();
    printf("\n");
}

void
LatticeNeumannCouplingNode :: computeLoadCouplingContribution(FloatArray &answer, TimeStep *stepN) {
    int nCouplingNodes = couplingNodes.giveSize();

    FloatArray contribution;
    answer.resize(this->directionVector.giveSize() );

    IntArray coupledModels;
    double waterPressureNew = 0., waterPressureOld = 0.;

    if ( domain->giveEngngModel()->giveMasterEngngModel() ) {
        ( static_cast< StaggeredProblem * >( domain->giveEngngModel()->giveMasterEngngModel() ) )->giveCoupledModels(coupledModels);
        if ( coupledModels.at(2) != 0 && !( stepN->giveNumber() <= domain->giveEngngModel()->giveNumberOfFirstStep() ) ) {
            for ( int k = 0; k < nCouplingNodes; k++ ) {
                Node *coupledNode = domain->giveEngngModel()->giveMasterEngngModel()->giveSlaveProblem(coupledModels.at(2) )->giveDomain(1)->giveNode(couplingNodes(k) );
                const auto &couplingCoords = coupledNode->giveCoordinates();
                TimeStep *previousStep = domain->giveEngngModel()->giveMasterEngngModel()->givePreviousStep();

                waterPressureNew = coupledNode->giveDofWithID(P_f)->giveUnknown(VM_Total, stepN);
                if ( !stepN->givePreviousStep()->isTheFirstStep() ) {
                    waterPressureOld = coupledNode->giveDofWithID(P_f)->giveUnknown(VM_Total, previousStep);
                } else {
                    waterPressureOld = 0.;
                }

                double waterPressure = waterPressureNew - waterPressureOld;
                double dist = distance(couplingCoords, coordinates);
                double factor = waterPressure * dist;
                contribution = this->directionVector;
                contribution.times(factor);
                answer.add(contribution);
            }
        }
    }

    return;
}
}
