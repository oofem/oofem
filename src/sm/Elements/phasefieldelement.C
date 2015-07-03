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

#include "../sm/Elements/phasefieldelement.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/Materials/structuralms.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "domain.h"
#include "mathfem.h"
#include "timestep.h"
#include <cstdio>


namespace oofem {

PhaseFieldElement :: PhaseFieldElement( int i, Domain *aDomain ) 
{  
    ///@todo will be set by the cross section later
    internalLength = 6.0;
    criticalEnergy = 1.0e0;
    relaxationTime = 1.0;
};

void
PhaseFieldElement :: computeLocationArrayOfDofIDs( const IntArray &dofIdArray, IntArray &answer )
{
    // Routine to extract compute the location array an element given an dofid array.
    answer.clear();
    NLStructuralElement *el = this->giveElement();
    int k = 0;
    for(int i = 1; i <= el->giveNumberOfDofManagers(); i++) {
        DofManager *dMan = el->giveDofManager( i );
        for(int j = 1; j <= dofIdArray.giveSize( ); j++) {

            if(dMan->hasDofID( (DofIDItem) dofIdArray.at( j ) )) {
                Dof *d = dMan->giveDofWithID( dofIdArray.at( j ) );
                answer.followedBy( k + d->giveNumber( ) );
            }
        }
        k += dMan->giveNumberOfDofs( );
    }
}

void
    ///@todo Remove these functions Jim, they work identically to the one in Element.
PhaseFieldElement :: computeDisplacementUnknowns(FloatArray &answer, ValueModeType valueMode, TimeStep *stepN)
{
    IntArray dofIdArray;
    this->giveDofManDofIDMask_u(dofIdArray);
    Element :: computeVectorOf(dofIdArray, valueMode, stepN, answer);
}

void
PhaseFieldElement :: computeDamageUnknowns(FloatArray &answer, ValueModeType valueMode, TimeStep *stepN)
{
    IntArray dofIdArray;
    this->giveDofManDofIDMask_d(dofIdArray);
    Element :: computeVectorOf(dofIdArray, valueMode, stepN, answer);
}

void
PhaseFieldElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // Computes the internal forces corresponding to the two fields u & d
    IntArray IdMask_u, IdMask_d;
    this->giveDofManDofIDMask_u( IdMask_u );
    this->giveDofManDofIDMask_d( IdMask_d );
    this->computeLocationArrayOfDofIDs( IdMask_u, loc_u );
    this->computeLocationArrayOfDofIDs( IdMask_d, loc_d );
    
    int ndofs = this->computeNumberOfDofs();
    answer.resize( ndofs);
    answer.zero();

    FloatArray answer_u(0);
    FloatArray answer_d(0);

    this->giveInternalForcesVector_u(answer_u, tStep, useUpdatedGpRecord);
    this->giveInternalForcesVector_d(answer_d, tStep, useUpdatedGpRecord);
    answer.assemble(answer_u, loc_u);
    answer.assemble(answer_d, loc_d);
}

void
PhaseFieldElement :: giveInternalForcesVector_u(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // computes int_V ( B_u^t * BSgima_u ) * dV
    FloatArray NStress, BStress, vGenStress, NS;
    FloatMatrix N, B;
    NLStructuralElement *el = this->giveElement( );

    answer.clear();
    for ( auto &gp: el->giveIntegrationRule(0) ) {
        double dV = el->computeVolumeAround(gp);
            
        // compute generalized stress measure
        el->computeBmatrixAt(gp, B);
        this->computeBStress_u(BStress, gp, tStep, useUpdatedGpRecord);
        answer.plusProduct(B, BStress, dV);
    }
    
}

void
PhaseFieldElement :: giveInternalForcesVector_d(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // Computes int_V ( N^t *Nstress_d  +  B^t * g_c*l*grad(d)  ) dV
    FloatArray NStress, BStress, NS, a_d, grad_d;
    FloatMatrix N, B;
    NLStructuralElement *el = this->giveElement( );
    this->computeDamageUnknowns( a_d, VM_Total, tStep );

    for ( auto &gp: el->giveIntegrationRule(0) ) {
        double dV = el->computeVolumeAround(gp);
            
        // compute generalized stress measures
        this->computeNd_matrixAt( *gp->giveCoordinates(), N);
        computeNStress_d(NStress, gp, tStep, useUpdatedGpRecord);
        NS.beTProductOf(N, NStress);
        answer.add(dV, NS);

        this->computeBd_matrixAt(gp, B );
        grad_d.beProductOf(B, a_d);
        double l = this->giveInternalLength();
        double g_c = this->giveCriticalEnergy( );
        BStress = grad_d * l * g_c;
        answer.plusProduct(B, BStress, dV);
      
    }
    
}

void
PhaseFieldElement :: computeBStress_u(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, int useUpdatedGpRecord)
{
    // computes G(d)*sig(u)
    NLStructuralElement *el = this->giveElement( );
    FloatArray strain, a_u;
    FloatMatrix B_u;

    el->computeBmatrixAt(gp, B_u);
    this->computeDisplacementUnknowns(a_u, VM_Total, tStep);
    strain.beProductOf(B_u, a_u);
    el->computeStressVector(answer, strain, gp, tStep);
    answer.times( this->computeG(gp, VM_Total, tStep) );

}

double
PhaseFieldElement :: computeFreeEnergy(GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
    FloatArray strain, stress;
    stress = matStat->giveTempStressVector();
    strain = matStat->giveTempStrainVector();
    return 0.5 * stress.dotProduct( strain );
}

void
PhaseFieldElement :: computeNStress_d(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, int useUpdatedGpRecord)
{
    // computes (t*/Delta t) *(d - d_old) + g_c/l*d + G'*Psibar

    //PhaseFieldCrossSection *cs = static_cast... 
    double Delta_t = tStep->giveTimeIncrement();
    double t_star = this->giveRelaxationTime();
    double d = computeDamageAt( gp, VM_Total, tStep ); 
    double Delta_d = computeDamageAt( gp, VM_Incremental, tStep );

    double l = this->giveInternalLength();
    double g_c = this->giveCriticalEnergy();
    double Gprim = this->computeGPrim(gp, VM_Total, tStep);
    double Psibar = this->computeFreeEnergy( gp, tStep );
    answer.resize( 1 );
    answer.at( 1 ) = t_star / Delta_t * Delta_d + g_c / l *d + Gprim * Psibar;
}

double 
PhaseFieldElement :: computeDamageAt(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // d = N_d * a_d
    NLStructuralElement *el = this->giveElement( );
    FloatArray dVec;
    computeDamageUnknowns(dVec, valueMode, stepN);
    FloatArray Nvec;
    el->giveInterpolation()->evalN(Nvec, *gp->giveCoordinates(), FEIElementGeometryWrapper(el));
    return Nvec.dotProduct(dVec);
}

double
PhaseFieldElement :: computeG(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // computes Dg/Dd = (1-d)^2 + r0
    double d = this->computeDamageAt(gp, valueMode, stepN);
    double r0 = 1.0e-10;
    return (1.0 - d) * (1.0 - d) + r0;
}

double 
PhaseFieldElement :: computeGPrim(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // compute -2*(1-d)
    double d = this->computeDamageAt(gp, valueMode, stepN);
    return -2.0 * (1.0 - d);
}

void
PhaseFieldElement :: computeNd_matrixAt(const FloatArray &lCoords, FloatMatrix &N)
{
    NLStructuralElement *el = this->giveElement( );
    FloatArray Nvec;
    el->giveInterpolation( )->evalN( Nvec, lCoords, FEIElementGeometryWrapper( el ) );
    N.resize(1, Nvec.giveSize());
    N.beNMatrixOf(Nvec,1);

}

void
PhaseFieldElement :: computeBd_matrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    // Returns the [numSpaceDim x nDofs] gradient matrix {B_d} of the receiver,
    // evaluated at gp.

    NLStructuralElement *el = this->giveElement( );
    FloatMatrix dNdx;
    el->giveInterpolation( )->evaldNdx( dNdx, *gp->giveCoordinates( ), FEIElementGeometryWrapper( el ) );
    answer.beTranspositionOf( dNdx );
}

int
PhaseFieldElement :: computeNumberOfDofs()
{
    ///@todo This is NOT correct procedure. This function may absolutely not ask the dof managers _anything_.
    NLStructuralElement *el = this->giveElement( );
    int nDofs = 0;
    for( int i = 1; i <= el->giveNumberOfDofManagers(); i++)
    {
        nDofs += el->giveDofManager( i )->giveNumberOfDofs();
    }
    return nDofs;
}

void
PhaseFieldElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
    ///@todo this part is enough to do once
    IntArray IdMask_u, IdMask_d;
    this->giveDofManDofIDMask_u( IdMask_u );
    this->giveDofManDofIDMask_d( IdMask_d );
    this->computeLocationArrayOfDofIDs( IdMask_u, loc_u );
    this->computeLocationArrayOfDofIDs( IdMask_d, loc_d );

    int nDofs = this->computeNumberOfDofs();
    answer.resize( nDofs, nDofs );
    answer.zero();

    FloatMatrix answer1, answer2, answer3, answer4;
    this->computeStiffnessMatrix_uu(answer1, rMode, tStep);
    this->computeStiffnessMatrix_ud(answer2, rMode, tStep);
    //this->computeStiffnessMatrix_du(answer3, rMode, tStep); //symmetric
    answer3.beTranspositionOf( answer2 );
    this->computeStiffnessMatrix_dd(answer4, rMode, tStep);
    
    answer.assemble( answer1, loc_u, loc_u );
    answer.assemble( answer2, loc_u, loc_d );
    answer.assemble( answer3, loc_d, loc_u );
    answer.assemble( answer4, loc_d, loc_d );
}

void
PhaseFieldElement :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    // This is the regular stiffness matrix times G
    FloatMatrix B, DB, N, D_B;
    NLStructuralElement *el = this->giveElement( );
    StructuralCrossSection *cs = dynamic_cast<StructuralCrossSection* > (el->giveCrossSection() );

    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);
    answer.clear();

    for ( auto gp: el->giveIntegrationRule(0) ) {
        double dV = el->computeVolumeAround(gp);
        // compute int_V ( B^t * D_B * B )dV
        el->computeBmatrixAt(gp, B );
        cs->giveCharMaterialStiffnessMatrix(D_B, rMode, gp, tStep);
        D_B.times( computeG(gp, VM_Total, tStep) );
        DB.beProductOf(D_B, B);

        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(B, DB, dV);
        } else {
            answer.plusProductUnsym(B, DB, dV);
        }    
        
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
    
}

void
PhaseFieldElement :: computeStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    FloatMatrix B, DB, N_d, DN, S(3,1);
    FloatArray stress;
    NLStructuralElement *el = this->giveElement( );
    StructuralCrossSection *cs = dynamic_cast<StructuralCrossSection* > (el->giveCrossSection() );

    answer.clear();

    for ( auto &gp: el->giveIntegrationRule(0) ) {
        StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );

        double dV = el->computeVolumeAround(gp);
        // compute int_V ( B^t * D_B * B )dV
        
        this->computeNd_matrixAt(*gp->giveCoordinates(), N_d);

        // stress   
        FloatArray reducedStrain, a_u;
        FloatMatrix B_u;
        el->computeBmatrixAt( gp, B_u );
        stress = matStat->giveTempStressVector();
        stress.times( this->computeGPrim( gp, VM_Total, tStep ) );
        
        S.setColumn(stress,1);
        DN.beProductOf(S, N_d);
        answer.plusProductUnsym(B_u, DN, dV);        
    }

}

void
PhaseFieldElement :: computeStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    double Delta_t = tStep->giveTimeIncrement();
    double t_star = this->giveRelaxationTime();
    double l = this->giveInternalLength();
    double g_c = this->giveCriticalEnergy();

    FloatMatrix B_d, N_d;
    //StructuralCrossSection *cs = dynamic_cast<StructuralCrossSection* > (this->giveCrossSection() );
    NLStructuralElement *el = this->giveElement( );
    answer.clear();

    for ( auto &gp: el->giveIntegrationRule(0) ) {
        double dV = el->computeVolumeAround(gp);
        
        this->computeNd_matrixAt(*gp->giveCoordinates(), N_d);
        this->computeBd_matrixAt(gp, B_d, 1, 3);
        
        double Gprim = this->computeGPrim(gp, VM_Total, tStep);
        
        double psiBar = this->computeFreeEnergy( gp, tStep );
        //double factorN = t_star/Delta_t + g_c/l + psiBar*Gbis;
        double factorN = t_star / Delta_t + g_c / l + psiBar*(-2.0);
        double factorB = g_c*l;
        
        answer.plusProductSymmUpper(N_d, N_d, factorN*dV);
        answer.plusProductSymmUpper(B_d, B_d, factorB*dV);   
    }

    answer.symmetrized();
}


IRResultType
PhaseFieldElement :: initializeFrom(InputRecord *ir)
{
    //IRResultType result;                // Required by IR_GIVE_FIELD macro
    //nlGeo = 0;

    return IRRT_OK;
}



//void
//PhaseFieldElement :: computeStressVectorAndLocalCumulatedStrain(FloatArray &answer, double localCumulatedStrain, GaussPoint *gp, TimeStep *stepN)
////PhaseFieldElement :: computeStressVector(FloatArray &answer, double localCumulatedStrain, GaussPoint *gp, TimeStep *stepN)
//{
//    NLStructuralElement *elem = this->giveNLStructuralElement();
//
//    nlGeo = elem->giveGeometryMode();
//    StructuralCrossSection *cs = this->giveNLStructuralElement()->giveStructuralCrossSection();
//
//    //this->computeDamage(alpha, gp, stepN);
//    if ( nlGeo == 0 ) {
//        FloatArray Epsilon;
//        this->computeLocalStrainVector(Epsilon, gp, stepN);
//        dpmat->giveRealStressVector(answer, gp, Epsilon, nlCumulatedStrain, stepN);
//    } else if ( nlGeo == 1 ) {
//        if ( elem->giveDomain()->giveEngngModel()->giveFormulation() == TL ) {
//            FloatArray vF;
//            this->computeDeformationGradientVector(vF, gp, stepN);
//            dpmat->giveFirstPKStressVector(answer, gp, vF, nlCumulatedStrain, stepN);
//        }
//
//    }
//}
//
//
//void
//PhaseFieldElement :: giveInternalForcesVector_u(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
//{
//    // computes int_V ( B_u^t * BSgima_u ) * dV
//    this->giveInternalForcesVectorGen(answer, tStep, useUpdatedGpRecord, 
//        NULL, 
//        &this->computeBmatrixAt,
//        NULL, 
//        &this->computeBStress_u,
//        &this->computeVolumeAround
//        );
//
//}
//
//void
//PhaseFieldElement :: giveInternalForcesVector_d(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
//{
//    // computes int_V ( N_d^t * NSgima_d + B_d^t * BSgima_d ) * dV
//    this->giveInternalForcesVectorGen(answer, tStep, useUpdatedGpRecord, 
//        &this->computeNmatrixAt, 
//        &this->computeBmatrixAt,
//        &this->computeNStress_d, 
//        &this->computeBStress_d,
//        &this->computeVolumeAround
//        );
//
//}
//
//void
//PhaseFieldElement :: computeBStress_d(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, int useUpdatedGpRecord)
//{
//    // computes g_c*l
//    // PhaseFieldCrossSection *cs = static_cast... 
//    double l = this->giveInternalLength();
//    double g_c = this->giveCriticalEnergy();
//    answer.resize(1);
//    answer.at(1) = g_c*l;
//}
//
//void
//PhaseFieldElement :: computeStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
//{
//
//}
//void
//PhaseFieldElement :: computeStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
//{
//   
//    StructuralCrossSection *cs = this->giveStructuralCrossSection();
//        
//    this->computeStiffnessMatrixGen(answer, rMode, tStep, 
//        NULL, 
//        &this->computeBmatrixAt,
//        NULL, 
//        &this->Duu_B,
//        this->computeVolumeAround 
//        );
//
//}
//
//void
//PhaseFieldElement :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
//{
//    // This is the regular stiffness matrix times G
//        
//    this->computeStiffnessMatrixGen(answer, rMode, tStep, 
//        NULL, 
//        &this->computeBmatrixAt,
//        NULL, 
//        &this->Duu_B,
//        this->computeVolumeAround 
//        );
//
//}
//
//void
//PhaseFieldElement :: Duu_B(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//{
//    // G*dSig/dEps
//    //StructuralCrossSection *cs = this->giveStructuralCrossSection();
//    //cs->giveCharMaterialStiffnessMatrix(answer, rMode, gp, tStep);
//
//    //answer.times( this->computeG(gp) );
//}





} // end namespace oofem
