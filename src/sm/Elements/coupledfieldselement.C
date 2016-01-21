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

#if 1
#include "../sm/Elements/coupledfieldselement.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/Elements/nlstructuralelement.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "cltypes.h"
#include "mathfem.h"
#include "unknownnumberingscheme.h"
#include <cstdio>

namespace oofem {
CoupledFieldsElement :: CoupledFieldsElement(int i, Domain *aDomain) : NLStructuralElement(i, aDomain)
// Constructor.
{
    nlGeo = 0;

}


void
CoupledFieldsElement :: computeLocationArrayOfDofIDs(const IntArray &dofIdArray, IntArray &answer)
{
    // Routine to extract compute the location array an element given an dofid array.
    answer.resize(0);
    
    int k = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        DofManager *dMan = this->giveDofManager(i);        
        for (int j = 1; j <= dofIdArray.giveSize(); j++ ) {   
            
            if ( dMan->hasDofID( (DofIDItem) dofIdArray.at(j) ) ) {
                Dof *d = dMan->giveDofWithID( dofIdArray.at(j) );
                answer.followedBy( k + d->giveNumber() );
                //answer.followedBy( k + j );
            }
        }
        k += dMan->giveNumberOfDofs( );
    }
}


void
CoupledFieldsElement :: computeVectorOfDofIDs(const IntArray &dofIdArray, ValueModeType valueMode, TimeStep *stepN, FloatArray &answer)
{
    // Routine to extract the solution vector for an element given an dofid array.
    // Size will be numberOfDofs and if a certain dofId does not exist a zero is used as value. 
    
    answer.resize( numberOfDofMans * dofIdArray.giveSize() ); // equal number of nodes for all fields
    answer.zero();
    int k = 1;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        DofManager *dMan = this->giveDofManager(i);        
        for (int j = 1; j <= dofIdArray.giveSize(); j++ ) {   
            
            if ( dMan->hasDofID( (DofIDItem) dofIdArray.at(j) ) ) {
                Dof *d = dMan->giveDofWithID( dofIdArray.at(j) );
                answer.at(k) = d->giveUnknown(valueMode, stepN);
            }
            k++;
        }
    }
}



void
CoupledFieldsElement :: giveInternalForcesVectorGen(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord, 
    void (*Nfunc)(GaussPoint*, FloatMatrix), void (*Bfunc)(GaussPoint*, FloatMatrix, int, int), //(GaussPoint*, FloatMatrix)
    void (*NStressFunc)(GaussPoint*, FloatArray), void (*BStressFunc)(GaussPoint*, FloatArray),
    double (*dVFunc)(GaussPoint*))
{
    // General implementation of internal forces that computes
    // f = sum_gp( N^T*GenStress_N + B^T*GenStress_B ) * dV

    FloatArray NStress, BStress, vGenStress, NS, BS;
    FloatMatrix N, B;

    for ( int j = 0; j < this->giveNumberOfIntegrationRules(); j++ ) {
        for ( auto &gp: this->giveIntegrationRule(j) ) {
            double dV  = this->computeVolumeAround(gp);
            
            // compute generalized stress measures
            if ( NStressFunc && Nfunc ) {
                Nfunc(gp, N);
                NStressFunc(gp, NStress);
                NS.beTProductOf(N, NStress);
                answer.add(dV, NS);
            }

            if ( BStressFunc && Bfunc ) {
                Bfunc(gp, B, 1, 3);
                BStressFunc(gp, BStress);
                BS.beTProductOf(B, BStress);
                answer.add(dV, BS);
            }

            
        }
    }
}




void
CoupledFieldsElement :: computeStiffnessMatrixGen(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep, 
        void (*Nfunc)(GaussPoint*, FloatMatrix), 
        void (*Bfunc)(GaussPoint*, FloatMatrix),
        void (*NStiffness)(FloatMatrix, MatResponseMode, GaussPoint*, TimeStep*), 
        void (*BStiffness)(FloatMatrix, MatResponseMode, GaussPoint*, TimeStep*),
        double (*volumeAround)(GaussPoint*) )
{
    FloatMatrix B, DB, N, DN, D_B, D_N;

    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answer.resize(0,0);

    for ( auto &gp: this->giveIntegrationRule(0) ) {
        double dV = this->computeVolumeAround(gp);


        // compute int_V ( N^t * D_N * N )dV
        if ( NStiffness && Nfunc ) {
            Nfunc(gp, N);
            NStiffness(D_N, rMode, gp, tStep);
            DN.beProductOf(D_N, N);
            if ( matStiffSymmFlag ) {
                answer.plusProductSymmUpper(N, DN, dV);
            } else {
                answer.plusProductUnsym(N, DN, dV);
            }
        }


        // compute int_V ( B^t * D_B * B )dV
        if ( BStiffness && Bfunc ) {
            Bfunc(gp, B);
            BStiffness(D_B, rMode, gp, tStep);
            DB.beProductOf(D_B, B);
            if ( matStiffSymmFlag ) {
                answer.plusProductSymmUpper(B, DB, dV);
            } else {
                answer.plusProductUnsym(B, DB, dV);
            }    
        }

    }


    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}



IRResultType
CoupledFieldsElement :: initializeFrom(InputRecord *ir)
{
    //IRResultType result;                // Required by IR_GIVE_FIELD macro
    //nlGeo = 0;

    return IRRT_OK;
}


} // end namespace oofem

#endif