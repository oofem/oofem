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

#include "structuralinterfaceelement.h"
#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"
#include "structuralinterfacecrosssection.h"
#include "feinterpol.h"
#include "domain.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "mathfem.h"



namespace oofem {
StructuralInterfaceElement :: StructuralInterfaceElement(int n, Domain *aDomain) : Element(n, aDomain)
{
    // Constructor. Creates an element with number n, belonging to aDomain.
    activityLtf = 0;
    initialDisplacements = NULL;
    nlGeometry = 0;
}


StructuralInterfaceElement :: ~StructuralInterfaceElement()
// Destructor.
{
    if ( initialDisplacements ) {
        delete initialDisplacements;
    }
}



void
StructuralInterfaceElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    // Computes the stiffness matrix of the receiver
    FloatMatrix N, D, DN;
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);

    answer.resize(0, 0);

    if ( !this->isActivated(tStep) ) {
        return;
    }

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    FloatMatrix rotationMatGtoL;
    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(j);

        if ( this->nlGeometry == 0 ) {
            this->giveStiffnessMatrix_Eng(D, rMode, ip, tStep);
        } else if ( this->nlGeometry == 1 ) {
            this->giveStiffnessMatrix_dTdj(D, rMode, ip, tStep);
        } else {
            OOFEM_ERROR("nlgeometry must be 0 or 1!")
        }

        this->computeTransformationMatrixAt(ip, rotationMatGtoL);
        D.rotatedWith(rotationMatGtoL, 't');                      // transform stiffness to global coord system

        this->computeNmatrixAt(ip, N);
        DN.beProductOf(D, N);
        double dA = this->computeAreaAround(ip);
        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(N, DN, dA);
        } else {
            answer.plusProductUnsym(N, DN, dA);
        }
    }


    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}


void
StructuralInterfaceElement :: computeSpatialJump(FloatArray &answer, IntegrationPoint *ip, TimeStep *tStep)
{
    // Computes the spatial jump vector (allways 3 components) at the Gauss point (ip) of
    // the receiver, at time step (tStep). jump = N*u
    FloatMatrix N;
    FloatArray u;

    if ( !this->isActivated(tStep) ) {
        answer.resize(3);
        answer.zero();
        return;
    }

    this->computeNmatrixAt(ip, N);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    answer.beProductOf(N, u);
}





void
StructuralInterfaceElement :: giveInternalForcesVector(FloatArray &answer,
                                                       TimeStep *tStep, int useUpdatedGpRecord)
{
    // Computes internal forces
    // if useGpRecord == 1 then data stored in ip->giveStressVector() are used
    // instead computing stressVector through this->ComputeStressVector();
    // this must be done after you want internal forces after element->updateYourself()
    // has been called for the same time step.

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    FloatMatrix N, rotationMatGtoL;
    FloatArray u, traction, tractionTemp, jump;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.resize(0);

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(i);
        this->computeNmatrixAt(ip, N);

        //if ( useUpdatedGpRecord == 1 ) {
        //    tractionTemp = static_cast<StructuralInterfaceMaterialStatus *> (ip->giveMaterialStatus())->giveTraction();
        //    // modify sizes if stored traction has 3 components and the spatial dimension is less
        //    int spatialDim = this->giveSpatialDimension();
        //    if ( traction.giveSize() == 3 && spatialDim == 1 ) {
        //        traction.setValues(1, tractionTemp.at(3) );
        //    } else if ( traction.giveSize() == 3 && spatialDim == 2 ) {
        //        traction.setValues(2, tractionTemp.at(1), tractionTemp.at(3) );
        //    } else {
        //        traction = tractionTemp;
        //    }
        //} else {
        if ( this->isActivated(tStep) ) {
            jump.beProductOf(N, u);
        } else {
            jump.resize(3);
            jump.zero();
        }
        this->computeTraction(traction, ip, jump, tStep);
        //}

        if ( traction.giveSize() == 0 ) {
            break;
        }

        // compute internal cohesive forces as f = N^T*traction dA
        double dA = this->computeAreaAround(ip);
        answer.plusProduct(N, traction, dA);
    }

    // if inactive update state, but no contribution to global system
    ///@todo why do we need to do this?
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}

void
StructuralInterfaceElement :: computeTraction(FloatArray &traction, IntegrationPoint *ip, FloatArray &jump, TimeStep *tStep)
{
    // Returns the traction in global coordinate system
    FloatMatrix rotationMatGtoL, F;
    this->computeTransformationMatrixAt(ip, rotationMatGtoL);
    jump.rotatedWith(rotationMatGtoL, 'n');      // transform jump to local coord system

    if ( this->nlGeometry == 0 ) {
        this->giveEngTraction(traction, ip, jump, tStep);
    } else if ( this->nlGeometry == 1 ) {
        ///@todo compute F in a proper way
        F.beUnitMatrix();
        F.rotatedWith(rotationMatGtoL, 'n');
        this->giveFirstPKTraction(traction, ip, jump, F, tStep);
    }

    traction.rotatedWith(rotationMatGtoL, 't');     // transform traction to global coord system
}



void
StructuralInterfaceElement :: giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep)
{
    // returns characteristics matrix of receiver according to mtrx

    if ( mtrx == StiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep);
    } else if ( mtrx == TangentStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep);
    } else if ( mtrx == SecantStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, SecantStiffness, tStep);
    } else if ( mtrx == ElasticStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, ElasticStiffness, tStep);
    } else {
        _error2( "giveCharacteristicMatrix: Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx) );
    }
}



void
StructuralInterfaceElement :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                                       TimeStep *tStep)
//
// returns characteristics vector of receiver according to mtrx
//
{
    if ( ( mtrx == InternalForcesVector ) && ( mode == VM_Total ) ) {
        this->giveInternalForcesVector(answer, tStep);
    } else if ( ( mtrx == LastEquilibratedInternalForcesVector ) && ( mode == VM_Total ) ) {
        /* here tstep is not relevant, we set useUpdatedGpRecord = 1
         * and this will cause to integrate internal forces using existing (nontemp, equlibrated) stresses in
         * statuses. Mainly used to compute reaction forces */
        this->giveInternalForcesVector(answer, tStep, 1);
    }
}

void
StructuralInterfaceElement :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);

    // record initial displacement if element not active
    if ( activityLtf && !isActivated(tStep) ) {
        if ( !initialDisplacements ) {
            initialDisplacements = new FloatArray();
        }

        this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, * initialDisplacements);
    }
}


void
StructuralInterfaceElement :: updateInternalState(TimeStep *tStep)
{
    // Updates the receiver at end of step.

    FloatArray tractionG, jumpL;
    FloatMatrix rotationMatGtoL;

    // force updating strains & stresses
    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        IntegrationRule *iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            IntegrationPoint *ip = iRule->getIntegrationPoint(j);
            this->computeSpatialJump(jumpL, ip, tStep);
            this->computeTraction(tractionG, ip, jumpL, tStep);
        }
    }
}

int
StructuralInterfaceElement :: checkConsistency()
{
    // check if the cross section is a 'StructuralInterfaceCrossSection'
    int result = 1;
    if ( !this->giveInterfaceCrossSection()->testCrossSectionExtension(this->giveInterfaceCrossSection()->crossSectionType) ) {
        OOFEM_ERROR2( "StructuralInterfaceElement :: checkConsistency : cross section %s is not a structural interface cross section", this->giveCrossSection()->giveClassName() );
        result = 0;
    }
    return result;
}


int
StructuralInterfaceElement :: giveIPValue(FloatArray &answer, IntegrationPoint *aIntegrationPoint, InternalStateType type, TimeStep *tStep)
{
    return Element :: giveIPValue(answer, aIntegrationPoint, type, tStep);
}





IRResultType
StructuralInterfaceElement :: initializeFrom(InputRecord *ir)
{
    return Element :: initializeFrom(ir);
}

void StructuralInterfaceElement :: giveInputRecord(DynamicInputRecord &input)
{
    Element :: giveInputRecord(input);
}


StructuralInterfaceCrossSection *StructuralInterfaceElement :: giveInterfaceCrossSection()
{
    return static_cast< StructuralInterfaceCrossSection * >( this->giveCrossSection() );
};



void
StructuralInterfaceElement :: giveEngTraction(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
{
    // Default implementation uses the First PK version if available.
    ///@todo Should probably convert to 2nd PK traction to be consistent with continuum elements.
    ///@todo compute F properly
    FloatMatrix F(3, 3);
    F.beUnitMatrix();
    this->giveFirstPKTraction(answer, gp, jump, F, tStep);
    // T_PK2 = inv(F)* T_PK1
}


void
StructuralInterfaceElement :: giveStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, IntegrationPoint *ip, TimeStep *tStep)
{
    // Default implementation uses the First PK version if available
    this->giveStiffnessMatrix_dTdj(answer, rMode, ip, tStep);
    ///@todo dT1dj = d(F*T2)/dj = dF/dj*T2 + F*dT2/dj - should we assume dFdj = 0? Otherwise it will be very complex!
}
} // end namespace oofem
