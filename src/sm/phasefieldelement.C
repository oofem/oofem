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
#if 0
#include "phasefieldelement.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "domain.h"
#include "mathfem.h"
#include "structuralcrosssection.h"
#include "timestep.h"
#include "structuralms.h"
#include <cstdio>



namespace oofem {

FEI2dQuadLin PhaseFieldElement :: interpolation_u(1, 2);
FEI2dQuadLin PhaseFieldElement :: interpolation_d(1, 2);

void
PhaseFieldElement :: computeDisplacementUnknowns(FloatArray &answer, ValueModeType valueMode, TimeStep *stepN)
{
    IntArray dofIdArray;
    this->giveDofManDofIDMask_u(dofIdArray);
    this->computeVectorOfDofIDs(dofIdArray, valueMode, stepN, answer);
}

void
PhaseFieldElement :: computeDamageUnknowns(FloatArray &answer, ValueModeType valueMode, TimeStep *stepN)
{
    IntArray dofIdArray;
    this->giveDofManDofIDMask_u(dofIdArray);
    this->computeVectorOfDofIDs(dofIdArray, valueMode, stepN, answer);
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
//        
//        if(cs->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
//            //dpmat->giveRealStressVector(answer, gp, Epsilon, nlCumulatedStrain, stepN);
//            return;
//        } else {
//            OOFEM_ERROR("giveRealStresses : unsupported mode");
//        }
//    } else if ( nlGeo == 1 ) {
//        if ( elem->giveDomain()->giveEngngModel()->giveFormulation() == TL ) {
//            FloatArray vF;
//            this->computeDeformationGradientVector(vF, gp, stepN);
//            if( cs->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
//              //  dpmat->giveFirstPKStressVector(answer, gp, vF, nlCumulatedStrain, stepN);
//                return;
//            } else {
//                OOFEM_ERROR("giveRealStresses : unsupported mode");
//            }
//        }
//
//    }
//}

void
PhaseFieldElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    //answer.resize(totalSize);
    answer.zero();
    FloatArray answer_u(0);
    FloatArray answer_d(0);

    this->giveInternalForcesVector_u(answer_u, tStep, useUpdatedGpRecord);
    this->giveInternalForcesVector_d(answer_d, tStep, useUpdatedGpRecord);
    IntArray IdMask_u, IdMask_d;
    this->giveDofManDofIDMask_u(IdMask_u);
    this->giveDofManDofIDMask_d(IdMask_d);
    this->computeLocationArrayOfDofIDs(IdMask_u, loc_u);
    this->computeLocationArrayOfDofIDs(IdMask_d, loc_d);
    answer.assemble(answer_u, loc_u);
    answer.assemble(answer_d, loc_d);
}

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
void
PhaseFieldElement :: giveInternalForcesVector_u(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // computes int_V ( B_u^t * BSgima_u ) * dV
     IntegrationRule *iRule = this->giveIntegrationRule(0);
    FloatArray NStress, BStress, vGenStress, NS, BS;
    FloatMatrix N, B;

    for ( int j = 0; j < this->giveNumberOfIntegrationRules(); j++ ) {
        IntegrationRule *iRule = this->giveIntegrationRule(j);

        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);
    
            double dV  = this->computeVolumeAround(gp);
            
            // compute generalized stress measure
                this->computeBmatrixAt(gp, B, 1, 3);
                this->computeBStress_u(BStress, gp, tStep, useUpdatedGpRecord);
                //BStressFunc(gp, BStress);
                BS.beTProductOf(B, BStress);
                answer.add(dV, BS);      
        }
    }
}



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
void
PhaseFieldElement :: giveInternalForcesVector_d(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
 IntegrationRule *iRule = this->giveIntegrationRule(0);
    FloatArray NStress, BStress, vGenStress, NS, BS;
    FloatMatrix N, B;

    for ( int j = 0; j < this->giveNumberOfIntegrationRules(); j++ ) {
        IntegrationRule *iRule = this->giveIntegrationRule(j);

        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);
    
            double dV  = this->computeVolumeAround(gp);
            
            // compute generalized stress measures
            //if ( NStressFunc && Nfunc ) {
                //Nfunc(gp, N);
                this->computeNd_matrixAt( *gp->giveCoordinates(), N);
                //NStressFunc(gp, NStress);
                computeNStress_d(NStress, gp, tStep, useUpdatedGpRecord);
                NStress.beTProductOf(N, NStress);
                answer.add(dV, NS);
            //}

            //if ( BStressFunc && Bfunc ) {
                //Bfunc(gp, B, 1, 3);
                this->computeBd_matrixAt(gp, B, 1, 3);
                //BStressFunc(gp, BStress);
                computeBStress_d(BStress, gp, tStep, useUpdatedGpRecord);
                BS.beTProductOf(B, BStress);
                answer.add(dV, BS);
            //}

            
        }
    }
}



void
PhaseFieldElement :: computeBStress_u(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, int useUpdatedGpRecord)
{
    // computes G(d)*sig(u)
    StructuralCrossSection *cs = dynamic_cast< StructuralCrossSection *> ( this->giveCrossSection() );
    FloatArray reducedStrain, a_u;
    FloatMatrix B_u;
    this->computeBmatrixAt(gp, B_u, 1, 3);

    this->computeDisplacementUnknowns(a_u, VM_Total, tStep);
    reducedStrain.beProductOf(B_u, a_u);
    cs->giveRealStresses(answer, gp, reducedStrain, tStep);

    answer.times( this->computeG(gp, VM_Total, tStep) );

}

void
PhaseFieldElement :: computeNStress_d(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, int useUpdatedGpRecord)
{
    // computes (t*/Delta t) *(d - d_old) + g_c/l*d + G'*Psibar
    
    //PhaseFieldCrossSection *cs = static_cast... 
    double Delta_t = tStep->giveTimeIncrement();
    // t_star = cs->give
    double t_star = 1.0;
    //double d = computeDamageAt(gp, tStep)
    //double d_old = computeDamageAt(gp, tStep_prev)
    double d = 0.0;
    double d_old = 0.0;
    // cs->giveInternalLength() 
    double l = 0.1;
    // cs->giveCriticalEnergy()
    double g_c = 1.0;
    // double Gprim = this->computeGPrim(gp);
    double Gprim = this->computeGPrim(gp, VM_Total, tStep);
    // double Psibar = cs->computeFreeEnergy(gp)
}

void
PhaseFieldElement :: computeBStress_d(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, int useUpdatedGpRecord)
{
    // computes g_c*l
    // PhaseFieldCrossSection *cs = static_cast... 
    double l = 0.1;
    double g_c = 1.0;
    answer.resize(1);
    answer.at(1) = g_c*l;
}

double 
PhaseFieldElement :: computeDamageAt(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // d = N_d * a_d
    FloatArray dVec;
    computeDamageUnknowns(dVec, valueMode, stepN);
    FloatArray Nvec;
    this->interpolation_d.evalN(Nvec, *gp->giveCoordinates(), FEIElementGeometryWrapper(this));
    return Nvec.dotProduct(dVec);
}

double
PhaseFieldElement :: computeG(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // compute (1-d)^2 + r0
    double d = this->computeDamageAt(gp, valueMode, stepN);
    double r0 = 1.0e-10;
    return (1.0 - d) * (1.0 - d) + r0;
}

double 
PhaseFieldElement :: computeGPrim(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // compute -2*d
    double d = this->computeDamageAt(gp, valueMode, stepN);
    return -2.0*d;
}

void
PhaseFieldElement :: computeNd_matrixAt(const FloatArray &lCoords, FloatMatrix &N)
{
    FloatArray Nvec;
    this->interpolation_d.evalN(Nvec, lCoords, FEIElementGeometryWrapper(this));
    N.resize(1, Nvec.giveSize());
    N.beNMatrixOf(Nvec,1);

}

void
PhaseFieldElement :: computeBd_matrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
//
// Returns the [2x4] gradient matrix {B_d} of the receiver,
// evaluated at aGaussPoint.
// a = ( d1, d2, d3, d4)
{

    FloatMatrix dnx;
    this->interpolation_d.evaldNdx( dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 4);
    answer.zero();

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, i) = dnx.at(i, 1);
        answer.at(2, i) = dnx.at(i, 2);
    }
}


void
PhaseFieldElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
    //this->setDisplacementLocationArray(locU, nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    //this->setNonlocalLocationArray(locK, nPrimNodes, nPrimVars, nSecNodes, nSecVars);


    //answer.resize(totalSize, totalSize);
    answer.zero();

    FloatMatrix answer1, answer2, answer3, answer4;
    this->computeStiffnessMatrix_uu(answer1, rMode, tStep);
    this->computeStiffnessMatrix_ud(answer2, rMode, tStep);
    this->computeStiffnessMatrix_du(answer3, rMode, tStep);
    this->computeStiffnessMatrix_dd(answer4, rMode, tStep);
    answer.assemble(answer1, loc_u);
    answer.assemble(answer2, loc_u, loc_d);
    answer.assemble(answer3, loc_d, loc_u);
    answer.assemble(answer4, loc_d);
}

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
void
PhaseFieldElement :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    // This is the regular stiffness matrix times G
    FloatMatrix B, DB, N, DN, D_B, D_N;
    StructuralCrossSection *cs = dynamic_cast<StructuralCrossSection* > (this->giveCrossSection() );

    IntegrationRule *iRule = this->giveIntegrationRule(0);
    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);
    answer.resize(0,0);

    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(j);
        
        double dV = this->computeVolumeAround(gp);
        // compute int_V ( B^t * D_B * B )dV
        this->computeBmatrixAt(gp, B, 1, 3);
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
    FloatMatrix B, DB, N, DN, D_B, D_N, S(3,1);
    FloatArray stress;
    StructuralCrossSection *cs = dynamic_cast<StructuralCrossSection* > (this->giveCrossSection() );

    IntegrationRule *iRule = this->giveIntegrationRule(0);
    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);
    answer.resize(0,0);

    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(j);
        
        double dV = this->computeVolumeAround(gp);
        // compute int_V ( B^t * D_B * B )dV
        this->computeBmatrixAt(gp, B, 1, 3);
        this->computeNd_matrixAt(*gp->giveCoordinates(), N);
        N.times( this->computeGPrim(gp, VM_Total, tStep) );
        stress = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        
        S.setColumn(stress,1);
        DN.beProductOf(S, N);
        answer.plusProductUnsym(B, DN, dV);        
    }

}

void
PhaseFieldElement :: computeStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

}

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
void
PhaseFieldElement :: computeStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

        double Delta_t = tStep->giveTimeIncrement();
    // t_star = cs->give
    double t_star = 1.0;
    //double d = computeDamageAt(gp, tStep)
    //double d_old = computeDamageAt(gp, tStep_prev)
    double d = 0.0;
    double d_old = 0.0;
    // cs->giveInternalLength() 
    double l = 0.1;
    // cs->giveCriticalEnergy()
    double g_c = 1.0;
    // double Gprim = this->computeGPrim(gp);
    

    FloatMatrix B, DB, N, DN, D_B, D_N, S(3,1);
    FloatArray stress;
    StructuralCrossSection *cs = dynamic_cast<StructuralCrossSection* > (this->giveCrossSection() );

    IntegrationRule *iRule = this->giveIntegrationRule(0);
    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);
    answer.resize(0,0);

    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(j);
        
        double dV = this->computeVolumeAround(gp);
        
        this->computeNd_matrixAt(*gp->giveCoordinates(), N);
        this->computeBd_matrixAt(gp, B, 1, 3);
        
        double Gprim = this->computeGPrim(gp, VM_Total, tStep);
        double freeEnergy = 1.0e-5;
        //double factorN = t_star/Delta_t + g_c/l + freeEnergy*Gprim;
        double factorN = t_star/Delta_t + g_c/l + freeEnergy*(-2.0);
        double factorB = l*l;
        
        answer.plusProductSymmUpper(N, N, factorN*dV);
        answer.plusProductSymmUpper(B, B, factorB*dV);   
          
    }

    answer.symmetrized();
}

void
PhaseFieldElement :: Duu_B(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    // G*dSig/dEps
    //StructuralCrossSection *cs = this->giveStructuralCrossSection();
    //cs->giveCharMaterialStiffnessMatrix(answer, rMode, gp, tStep);

    //answer.times( this->computeG(gp) );
}


IRResultType
PhaseFieldElement :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                // Required by IR_GIVE_FIELD macro
    //nlGeo = 0;

    return IRRT_OK;
}



} // end namespace oofem

#endif