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

#include "latticedirichletcouplingnode.h"
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

  REGISTER_DofManager(LatticeDirichletCouplingNode);

LatticeDirichletCouplingNode :: LatticeDirichletCouplingNode(int n, Domain *aDomain) :
    Node(n, aDomain)
{}

LatticeDirichletCouplingNode :: ~LatticeDirichletCouplingNode()
// Destructor.
{}

void
LatticeDirichletCouplingNode :: initializeFrom(InputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
    FloatArray triplets;

    Node :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, couplingElements, _IFT_LatticeDirichletCouplingNode_couplingelements);

}

void LatticeDirichletCouplingNode :: printYourself()
// Prints the receiver on screen.
{
    int i;
    double x, y;

    x = this->giveCoordinate(1);
    y = this->giveCoordinate(2);
    printf("LatticeDirichletCouplingNode %d    coord : x %f  y %f\n", number, x, y);
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
LatticeDirichletCouplingNode :: giveUnknownVector(FloatArray &answer, const IntArray &dofIDArry, ValueModeType mode, TimeStep *tStep, bool padding)
//giveUnknownVector(FloatArray &answer, const IntArray &dofMask, ValueModeType mode, TimeStep *stepN)
{

    int size;
    IntArray dofArray;

    answer.resize( size = dofIDArry.giveSize() );

    for ( int i = 1; i <= size; i++ ) {
        auto pos = this->findDofWithDofId( ( DofIDItem ) dofIDArry.at(i) );
#ifdef DEBUG
        if ( pos == this->end() ) {
            OOFEM_ERROR("Couldn't find dof with Dof ID %d", dofIDArry.at(i));
        }
#endif
        answer.at(i) = (*pos)->giveUnknown(mode, tStep);
	activedirichletbc = computeUnknownCouplingContribution(tStep);
	answer.at(i) += activedirichletbc;

    }

    // Transform to global c.s.
    FloatMatrix L2G;
    if ( this->computeL2GTransformation(L2G, dofIDArry) ) {
        answer.rotatedWith(L2G, 'n');
    }


    // if ( !hasSlaveDofs ) {
    //     int j, size;
    //     IntArray dofArray;
    //     this->activedirichletbc = 0;
    //     answer.resize( size = dofMask.giveSize() );
    //     this->giveDofArray(dofMask, dofArray);

    //     for ( j = 1; j <= size; j++ ) {
    //         answer.at(j) = this->giveDof( dofArray.at(j) )->giveUnknown(mode, stepN);

    //         activedirichletbc = computeUnknownCouplingContribution(stepN);

    //         answer.at(j) += activedirichletbc;
    //     }
    // } else {
    //     int i, k, indx;
    //     IntArray dofArray;
    //     FloatArray mstrUnknwns;

    //     this->giveDofArray(dofMask, dofArray);
    //     answer.resize( giveNumberOfPrimaryMasterDofs(dofArray) );

    //     for ( k = 1, i = 1; i <= dofArray.giveSize(); i++ ) {
    //         indx = dofArray.at(i);
    //         if ( !this->giveDof(indx)->isPrimaryDof() ) { // slave DOF
    //             this->giveDof(indx)->giveUnknowns(mstrUnknwns, mode, stepN);
    //             answer.copySubVector(mstrUnknwns, k);
    //             k += mstrUnknwns.giveSize();
    //         } else {
    //             answer.at(k++) = this->giveDof(indx)->giveUnknown(mode, stepN);
    //         }
    //     }
    // }
}

double
LatticeDirichletCouplingNode :: computeUnknownCouplingContribution(TimeStep *stepN) {
    double distance = 0.;
    double result = 0.;

    int nCouplingElements = couplingElements.giveSize();

    FloatArray couplingGpCoords;
    double normalStress;

    double nominator = 0;
    double denominator = 0;

#ifdef __SM_MODULE
    IntArray coupledModels;
    if ( domain->giveEngngModel()->giveMasterEngngModel() ) {
        ( static_cast< StaggeredProblem * >( domain->giveEngngModel()->giveMasterEngngModel() ) )->giveCoupledModels(coupledModels);
        if ( coupledModels.at(1) != 0 && !(stepN->giveNumber()<=domain->giveEngngModel()->giveNumberOfFirstStep()) ) {
            LatticeStructuralElement *coupledElement;
            for ( int i = 1; i <= nCouplingElements; i++ ) {
                couplingGpCoords.zero();
                coupledElement  = static_cast< LatticeStructuralElement * >( domain->giveEngngModel()->giveMasterEngngModel()->giveSlaveProblem( coupledModels.at(1) )->giveDomain(1)->giveElement( couplingElements.at(i) ) );

                coupledElement->giveGpCoordinates(couplingGpCoords);
		if(coupledElement->hasBeenUpdated() == 1){
		  normalStress = coupledElement->giveOldNormalStress();
		}
		else{
		  normalStress = coupledElement->giveNormalStress();
		}
		  
		if(normalStress > 0.){
		  normalStress = 0.;
		}

                distance = sqrt( pow(couplingGpCoords.at(1) - coordinates.at(1), 2.) + pow(couplingGpCoords.at(2) - coordinates.at(2), 2.) );

                nominator = nominator + distance * normalStress;
                denominator = denominator + distance;
            }

            result = nominator / denominator;
        }

#endif
}

return result;
}

void LatticeDirichletCouplingNode :: printOutputAt(FILE *stream, TimeStep *stepN)
{
    int i;

#if defined( __PARALLEL_MODE ) || defined( __ENABLE_COMPONENT_LABELS )
    fprintf( stream, "%-8s%8d (%8d):\n", this->giveClassName(), this->giveLabel(), this->giveNumber() );
#else
    fprintf( stream, "%-8s%8d:\n", this->giveClassName(), this->giveGlobalNumber() );
#endif

    for ( i = 1; i <= this->giveNumberOfDofs(); i++ ) {
        fprintf(stream, "  dof %d f % .16e\n", this->giveGlobalNumber(), activedirichletbc);
    }
}
  
} // end namespace oofem
