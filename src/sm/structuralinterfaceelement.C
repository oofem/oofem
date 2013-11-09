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

#if 0
void
StructuralInterfaceElement :: computeConstitutiveMatrixAt(FloatMatrix &answer,
                                                 MatResponseMode rMode, IntegrationPoint *ip,
                                                 TimeStep *tStep)
// Returns the  material matrix {E= dT/dJ} of the receiver.
// rMode parameter determines type of stiffness matrix to be requested
// (tangent, secant, ...)
{
    if ( this->nlGeometry == 0 ) {
        this->giveInterfaceCrossSection()->give3dStiffnessMatrix_Eng(answer, rMode, ip, tStep);
    } else if ( this->nlGeometry = 1 ) {
        this->giveInterfaceCrossSection()->give3dStiffnessMatrix_dTdj(answer, rMode, ip, tStep);
    } else {

        OOFEM_ERROR("Spatial dimension must be 1-3!")
    }

}
#endif




void
StructuralInterfaceElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    // Computes the stiffness matrix of the receiver
    FloatMatrix N, D, DN;
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode, this->material);

    answer.resize(0, 0);

    if ( !this->isActivated(tStep) ) {
        return;
    }

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    FloatMatrix rotationMatGtoL;
    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(j);
        
        if ( this->nlGeometry == 0 ) {
            this->giveInterfaceCrossSection()->give3dStiffnessMatrix_Eng(D, rMode, ip, tStep);
        } else if ( this->nlGeometry = 1 ) {
            this->giveInterfaceCrossSection()->give3dStiffnessMatrix_dTdj(D, rMode, ip, tStep);
        } else {
            OOFEM_ERROR("nlgemoetry must be 0 or 1!")
        }
        //this->computeConstitutiveMatrixAt(D, rMode, ip, tStep);

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
        u.subtract(*initialDisplacements);
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
    StructuralInterfaceCrossSection *CS = this->giveInterfaceCrossSection();

    FloatMatrix N, rotationMatGtoL;
    FloatArray u, traction, jump;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(*initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.resize(0);

    FloatMatrix F(3,3);
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(i);
        this->computeNmatrixAt(ip, N);

        if ( useUpdatedGpRecord == 1 ) {
            traction = this->giveInterfaceCrossSection()->giveTraction(ip);
        } else {
            if ( this->isActivated(tStep) ) {
                jump.beProductOf(N, u);
            } else {
                jump.resize(3);
                jump.zero();
            }
            this->computeTransformationMatrixAt(ip, rotationMatGtoL); 
            //jump.printYourself();
            jump.rotatedWith(rotationMatGtoL, 'n');             // transform jump to local coord system
            //jump.printYourself();
            if( this->nlGeometry == 0 ) {
                CS->giveEngTraction_3d(traction, ip, jump, tStep); 
            } else if( this->nlGeometry == 1) {
                //compute F
                F.beUnitMatrix();
                F.rotatedWith(rotationMatGtoL, 'n');
                CS->giveFirstPKTraction_3d(traction, ip, jump, F, tStep); 
            }

            traction.rotatedWith(rotationMatGtoL,'t');          // transform traction to global coord system
        }

        if ( traction.giveSize() == 0 ) {
            break;
        }

        // compute internal cohesive forces as f = N^T*traction dA
        double dA = this->computeAreaAround(ip);
        answer.plusProduct(N, traction, dA);
    }
    answer.printYourself();
    // if inactive update state, but no contribution to global system
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}

void 
StructuralInterfaceElement :: computeTraction(FloatArray &traction, IntegrationPoint *ip, FloatArray &jump, TimeStep *tStep)
{
    StructuralInterfaceCrossSection *CS = this->giveInterfaceCrossSection();
/*    if ( this->giveSpatialDimension() == 1 ) {
        CS->giveEngTraction_1d(traction, ip, jump,  tStep);
    } else if ( this->giveSpatialDimension() == 2 ) {
        CS->giveEngTraction_2d(traction, ip, jump,  tStep);
    } else if ( this->giveSpatialDimension() == 3 ) {
    */
    CS->giveEngTraction_3d(traction, ip, jump,  tStep);
    //} else {
        // Will it ever run this part?
    //    OOFEM_ERROR("Spatial dimension must be 1-3!")
    //}
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

        this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, *initialDisplacements);
    }
}


void
StructuralInterfaceElement :: updateInternalState(TimeStep *tStep)
{
    // Updates the receiver at end of step.

    FloatArray traction, jump;
    FloatMatrix rotationMatGtoL;

    // force updating strains & stresses
    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        IntegrationRule *iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            IntegrationPoint *ip = iRule->getIntegrationPoint(j);
            this->computeSpatialJump(jump, ip, tStep);
            this->computeTransformationMatrixAt(ip, rotationMatGtoL); 
            jump.rotatedWith(rotationMatGtoL, 'n');             // transform jump to local coord system
            this->computeTraction(traction, ip, jump, tStep);
        }
    }
}

int
StructuralInterfaceElement :: checkConsistency()
{
    // check if the cross section is a 'StructuralInterfaceCrossSection'
    int result = 1;
    if ( !this->giveInterfaceCrossSection()->testCrossSectionExtension( this->giveInterfaceCrossSection()->crossSectionType ) ) {
        OOFEM_ERROR2("StructuralInterfaceElement :: checkConsistency : cross section %s is not a structural interface cross section", this->giveCrossSection()->giveClassName() );
        result = 0;
    }
    return result;
}


int
StructuralInterfaceElement :: giveIPValue(FloatArray &answer, IntegrationPoint *aIntegrationPoint, InternalStateType type, TimeStep *atTime)
{
    return Element :: giveIPValue(answer, aIntegrationPoint, type, atTime);
}





IRResultType
StructuralInterfaceElement :: initializeFrom(InputRecord *ir)
{
    return Element :: initializeFrom(ir);
}

void StructuralInterfaceElement :: giveInputRecord(DynamicInputRecord &input)
{
	Element :: giveInputRecord(input);

	/// TODO: Should initialDisplacements be stored? /ES
}


StructuralInterfaceCrossSection *StructuralInterfaceElement :: giveInterfaceCrossSection() 
{ 
    return static_cast< StructuralInterfaceCrossSection * > ( this->giveCrossSection() ); 
};

} // end namespace oofem
