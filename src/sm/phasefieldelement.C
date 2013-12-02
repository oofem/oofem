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
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "cltypes.h"
#include "structuralms.h"
#include "mathfem.h"
#include "structuralcrosssection.h"
#include "nlstructuralelement.h"
#include "nonlocalbarrier.h"
#include "nlstructuralelement.h"
#include "graddpmaterialextensioninterface.h"

#include <cstdio>

namespace oofem {
PhaseFieldElement :: PhaseFieldElement(int i, Domain *aDomain) : NLStructuralElement(i, aDomain)
// Constructor.
{
    nlGeo = 0;
    averType = 0;
}

void
PhaseFieldElement :: setDisplacementLocationArray(IntArray &answer, int nPrimNodes, int nPrimVars, int nSecNodes, int nSecVars)
{
    answer.resize(locSize);

    for ( int i = 1; i <= totalSize; i++ ) {
        if ( i < nSecNodes * nPrimVars + 1 ) {
            answer.at(i) = i + ( int )( ( ( i - 1 ) / nPrimVars ) ) * nSecVars;
        } else if ( i > nSecNodes * ( nPrimVars + nSecVars ) )  {
            answer.at(i - nSecVars * nSecNodes) = i;
        }
    }
}

void
PhaseFieldElement :: setNonlocalLocationArray(IntArray &answer, int nPrimNodes, int nPrimVars, int nSecNodes, int nSecVars)
{
    answer.resize(nlSize);
    for ( int i = 1; i <= nlSize; i++ ) {
        answer.at(i) = i * nPrimVars + i;
    }
}


void
PhaseFieldElement :: computeDisplacementDegreesOfFreedom(FloatArray &answer, TimeStep *stepN)
{
//    StructuralElement *elem = this->giveStructuralElement();
    FloatArray u;
    answer.resize(locSize);

  //  elem->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);
    u.resizeWithValues(totalSize);
    for ( int i = 1; i <= locSize; i++ ) {
        answer.at(i) = u.at( locU.at(i) );
    }
}

void PhaseFieldElement :: computeNonlocalDegreesOfFreedom(FloatArray &answer, TimeStep *stepN)
{
    //StructuralElement *elem = this->giveStructuralElement();
    FloatArray u;
    answer.resize(nlSize);

//    elem->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);
    u.resizeWithValues(totalSize);
    for  ( int i = 1; i <= nlSize; i++ ) {
        answer.at(i) = u.at( locK.at(i) );
    }
}

void
PhaseFieldElement :: computeStressVectorAndLocalCumulatedStrain(FloatArray &answer, double localCumulatedStrain, GaussPoint *gp, TimeStep *stepN)
//PhaseFieldElement :: computeStressVector(FloatArray &answer, double localCumulatedStrain, GaussPoint *gp, TimeStep *stepN)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();

    nlGeo = elem->giveGeometryMode();
    StructuralCrossSection *cs = this->giveNLStructuralElement()->giveStructuralCrossSection();
    GradDpMaterialExtensionInterface *dpmat = static_cast< GradDpMaterialExtensionInterface * >( 
            cs->giveMaterialInterface( GradDpMaterialExtensionInterfaceType, gp ) );

    if ( !dpmat ) {
        OOFEM_ERROR("PhaseFieldElement :: computeStiffnessMatrix_uu - Material doesn't implement the required DpGrad interface!");
    }

    //this->computeDamage(alpha, gp, stepN);
    if ( nlGeo == 0 ) {
        FloatArray Epsilon;
        this->computeLocalStrainVector(Epsilon, gp, stepN);
        
        if(cs->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
            //dpmat->giveRealStressVector(answer, gp, Epsilon, nlCumulatedStrain, stepN);
            return;
        } else {
            OOFEM_ERROR("giveRealStresses : unsupported mode");
        }
    } else if ( nlGeo == 1 ) {
        if ( elem->giveDomain()->giveEngngModel()->giveFormulation() == TL ) {
            FloatArray vF;
            this->computeDeformationGradientVector(vF, gp, stepN);
            if( cs->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
              //  dpmat->giveFirstPKStressVector(answer, gp, vF, nlCumulatedStrain, stepN);
                return;
            } else {
                OOFEM_ERROR("giveRealStresses : unsupported mode");
            }
        }

    }
}





void
PhaseFieldElement :: giveNonlocalInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    double dV, localCumulatedStrain=0.;
    GaussPoint *gp;
    NLStructuralElement *elem = this->giveNLStructuralElement();
    IntegrationRule *iRule =  elem->giveIntegrationRule(0);
    FloatMatrix stiffKappa, Nk;
    FloatArray fKappa(nlSize), aux(nlSize), dKappa, stress;

    aux.zero();
    int size = nSecVars * nSecNodes;

    //set displacement and nonlocal location array
    this->setDisplacementLocationArray(locU, nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    this->setNonlocalLocationArray(locK, nPrimNodes, nPrimVars, nSecNodes, nSecVars);

    answer.resize(size);
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        this->computeNkappaMatrixAt(gp, Nk);
        for ( int j = 1; j <= nlSize; j++ ) {
            fKappa.at(j) = Nk.at(1, j);
        }

        dV  = elem->computeVolumeAround(gp);
        this->computeStressVectorAndLocalCumulatedStrain(stress, localCumulatedStrain, gp, tStep);
        fKappa.times(localCumulatedStrain);
        fKappa.times(-dV);
        aux.add(fKappa);
    }

    this->computeStiffnessMatrix_kk(stiffKappa, TangentStiffness, tStep);
    this->computeNonlocalDegreesOfFreedom(dKappa, tStep);
    answer.beProductOf(stiffKappa, dKappa);
    answer.add(aux);
}


void
PhaseFieldElement :: giveInternalForcesVectorGen(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord, 
    void (*Nfunc)(GaussPoint*, FloatMatrix), void (*Bfunc)(GaussPoint*, FloatMatrix),
    void (*NStressFunc)(GaussPoint*, FloatArray), void (*BStressFunc)(GaussPoint*, FloatArray),
    double (*dVFunc)(GaussPoint*))
{
    
    IntegrationRule *iRule = this->giveIntegrationRule(0);
    FloatArray NStress, BStress, vGenStress, NS, BS;
    FloatMatrix N, B;

    for ( int j = 0; j < this->giveNumberOfIntegrationRules(); j++ ) {
        IntegrationRule *iRule = this->giveIntegrationRule(j);

        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);

            // Compute N and B matrices
            Nfunc(gp, N);
            Bfunc(gp, B);
            //this->computeBmatrixAt(gp, B);
            //this->computeNmatrixAt(gp, N);
    
            // compute generalized stress measures
            NStressFunc(gp, NStress);
            BStressFunc(gp, BStress);

            // Compute (generalized) internal forces as f = sum_gp( N^T*GenStress_N + B^T*GenStress_B ) * dV
            double dV  = this->computeVolumeAround(gp);
            NS.beTProductOf(N, NStress);
            BS.beTProductOf(B, BStress);
            answer.add(dV, NS);
            answer.add(dV, BS);
        }
    }
}


void
PhaseFieldElement :: giveLocalInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{

   // example -> (answer, tstep, upgp, N, B, sigN, sigB, dV)
    this->giveInternalForcesVectorGen(answer, tStep, useUpdatedGpRecord, 
        &this->computeNkappaMatrixAt, &this->computeBkappaMatrixAt,
        &this->computeNStressAt, &this->computeNStressAt, &this->computeVolumeAround);

}



void
PhaseFieldElement :: computeStiffnessMatrixGen(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    FloatMatrix B, D, DB;

    StructuralCrossSection *cs = this->giveNLStructuralElement()->giveStructuralCrossSection();
    IntegrationRule *iRule = elem->giveIntegrationRule(0);
    bool matStiffSymmFlag = elem->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answer.resize(locSize, locSize);
    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(j);
        
            // Compute N and B matrices
            //Nfunc(gp, N);
            //Bfunc(gp, B);
            

        //dpmat->givePDGradMatrix_uu(D, rMode, gp, tStep);
        double dV = elem->computeVolumeAround(gp);
        DB.beProductOf(D, B);
        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(B, DB, dV);
        } else {
            answer.plusProductUnsym(B, DB, dV);
        }
    }

    if ( elem->domain->giveEngngModel()->giveFormulation() == AL ) {
        FloatMatrix initialStressMatrix;
        elem->computeInitialStressMatrix(initialStressMatrix, tStep);
        answer.add(initialStressMatrix);
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}






void
PhaseFieldElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    answer.resize(totalSize);
    answer.zero();
    FloatArray answerU;
    answerU.resize(locSize);
    answer.zero();
    FloatArray answerK(nlSize);
    answerK.zero();

    this->giveNonlocalInternalForcesVector(answerK, tStep, useUpdatedGpRecord);
    this->giveLocalInternalForcesVector(answerU, tStep, useUpdatedGpRecord);


    answer.assemble(answerU, locU);
    answer.assemble(answerK, locK);
}

void
PhaseFieldElement :: computeForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
{
    //set displacement and nonlocal location array
    this->setDisplacementLocationArray(locU, nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    this->setNonlocalLocationArray(locK, nPrimNodes, nPrimVars, nSecNodes, nSecVars);


    FloatArray localForces(locSize);
    FloatArray nlForces(nlSize);
    answer.resize(totalSize);

    this->computeLocForceLoadVector(localForces, stepN, mode);

    answer.assemble(localForces, locU);
    answer.assemble(nlForces, locK);
}

void
PhaseFieldElement :: computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
{
    //set displacement and nonlocal location array
    this->setDisplacementLocationArray(locU, nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    this->setNonlocalLocationArray(locK, nPrimNodes, nPrimVars, nSecNodes, nSecVars);

    FloatArray localForces(locSize);
    FloatArray nlForces(nlSize);
    answer.resize(totalSize);
    answer.zero();

    this->computeLocNonForceLoadVector(answer, stepN, mode);
}


/************************************************************************/
void
PhaseFieldElement :: computeLocForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further sobstract part corresponding to non-nodeal loading.
{
    FloatMatrix T;
    NLStructuralElement *elem = this->giveNLStructuralElement();
    elem->computeLocalForceLoadVector(answer, stepN, mode);

    // transform result from global cs to nodal cs. if necessary
    if ( answer.isNotEmpty() ) {
        if ( elem->computeGtoLRotationMatrix(T) ) {
            // first back to global cs from element local
            answer.rotatedWith(T, 't');
        }
    } else   {
        answer.resize(locSize);
        answer.zero();
    }
}


void
PhaseFieldElement :: computeLocNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Computes the load vector of the receiver, at stepN.
{
    FloatArray helpLoadVector;
    StructuralElement *elem = this->giveStructuralElement();
    answer.resize(0);

    elem->computePrescribedStrainLoadVectorAt(helpLoadVector, stepN, mode);
    if ( helpLoadVector.giveSize() ) {
        answer.add(helpLoadVector);
    }
}


void
PhaseFieldElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
    this->setDisplacementLocationArray(locU, nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    this->setNonlocalLocationArray(locK, nPrimNodes, nPrimVars, nSecNodes, nSecVars);


    answer.resize(totalSize, totalSize);
    answer.zero();

    FloatMatrix answer1, answer2, answer3, answer4;
    this->computeStiffnessMatrix_uu(answer1, rMode, tStep);
    this->computeStiffnessMatrix_uk(answer2, rMode, tStep);
    this->computeStiffnessMatrix_ku(answer3, rMode, tStep);
    this->computeStiffnessMatrix_kk(answer4, rMode, tStep);
    answer.assemble(answer1, locU);
    answer.assemble(answer2, locU, locK);
    answer.assemble(answer3, locK, locU);
    answer.assemble(answer4, locK);
}




void
PhaseFieldElement :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();
    FloatMatrix B, D, DB;

    StructuralCrossSection *cs = this->giveNLStructuralElement()->giveStructuralCrossSection();

    nlGeo = elem->giveGeometryMode();
    IntegrationRule *iRule = elem->giveIntegrationRule(0);
    bool matStiffSymmFlag = elem->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answer.resize(locSize, locSize);
    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(j);
        
        GradDpMaterialExtensionInterface *dpmat = dynamic_cast< GradDpMaterialExtensionInterface * >( 
                    cs->giveMaterialInterface( GradDpMaterialExtensionInterfaceType, gp ) );
        if ( !dpmat ) {
            OOFEM_ERROR("PhaseFieldElement :: computeStiffnessMatrix_uk - Material doesn't implement the required DpGrad interface!");
        }
        if ( nlGeo == 0 ) {
            elem->computeBmatrixAt(gp, B);
        } else if ( nlGeo == 1 ) {
            if ( elem->domain->giveEngngModel()->giveFormulation() == AL ) {
                elem->computeBmatrixAt(gp, B);
            } else {
                elem->computeBHmatrixAt(gp, B);
            }
        }

        dpmat->givePDGradMatrix_uu(D, rMode, gp, tStep);
        double dV = elem->computeVolumeAround(gp);
        DB.beProductOf(D, B);
        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(B, DB, dV);
        } else {
            answer.plusProductUnsym(B, DB, dV);
        }
    }

    if ( elem->domain->giveEngngModel()->giveFormulation() == AL ) {
        FloatMatrix initialStressMatrix;
        elem->computeInitialStressMatrix(initialStressMatrix, tStep);
        answer.add(initialStressMatrix);
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}






void
PhaseFieldElement :: computeStiffnessMatrix_ku(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    double dV;
    NLStructuralElement *elem = this->giveNLStructuralElement();
    FloatMatrix B, DB, D, Nk, NkT, NkDB;
    StructuralCrossSection *cs = this->giveNLStructuralElement()->giveStructuralCrossSection();

    answer.resize(nlSize, locSize);

    IntegrationRule *iRule = elem->giveIntegrationRule(0);
    int nlGeo = this->giveNLStructuralElement()->giveGeometryMode();

    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(j);

        GradDpMaterialExtensionInterface *dpmat = dynamic_cast< GradDpMaterialExtensionInterface * >( 
                    cs->giveMaterialInterface( GradDpMaterialExtensionInterfaceType, gp ) );
        if ( !dpmat ) {
            OOFEM_ERROR("PhaseFieldElement :: computeStiffnessMatrix_uk - Material doesn't implement the required DpGrad interface!");
        }

        elem->computeBmatrixAt(gp, B);
        if ( nlGeo == 1 ) {
            if ( elem->domain->giveEngngModel()->giveFormulation() == AL ) {
                elem->computeBmatrixAt(gp, B);
            } else {
                elem->computeBHmatrixAt(gp, B);
            }
        }

        dpmat->givePDGradMatrix_ku(D, rMode, gp, tStep);
        this->computeNkappaMatrixAt(gp, Nk);
        NkT.beTranspositionOf(Nk);
        dV = elem->computeVolumeAround(gp);
        DB.beProductOf(D, B);
        NkDB.beProductOf(NkT, DB);
        NkDB.times(-dV);
        answer.add(NkDB);

        if ( dpmat->giveAveragingType() == 2 ) {
            double dl1, dl2, dl3;
            FloatArray Gk;
            FloatMatrix DB, LDB;
            FloatMatrix Bk, BkT, Nk, BktM22, BktM22Gk, BktM12, BktM12Gk, M22(2, 2), M12(2, 2);
            FloatMatrix dL1(1, 3), dL2(1, 3), result, result1, result2, dLdS, n(2, 2);

            this->computeBkappaMatrixAt(gp, Bk);
            BkT.beTranspositionOf(Bk);
            dpmat->givePDGradMatrix_uu(D, rMode, gp, tStep);
            dpmat->givePDGradMatrix_LD(dLdS, rMode, gp, tStep);
            this->computeNonlocalGradient(Gk, gp, tStep);

            dl1 = dLdS.at(3, 3);
            dl2 = dLdS.at(4, 4);
            dl3 = dLdS.at(5, 5);

            n.at(1, 1) = dLdS.at(1, 1);
            n.at(1, 2) = dLdS.at(1, 2);
            n.at(2, 1) = dLdS.at(2, 1);
            n.at(2, 2) = dLdS.at(2, 2);
            // first term BkT M22 G L1 D B
            // M22 = n2 \otimes n2
            M22.at(1, 1) = n.at(1, 2) * n.at(1, 2);
            M22.at(1, 2) = n.at(1, 2) * n.at(2, 2);
            M22.at(2, 1) = n.at(2, 2) * n.at(1, 2);
            M22.at(2, 2) = n.at(2, 2) * n.at(2, 2);
            // dL1
            dL1.at(1, 1) = dl1 * n.at(1, 1) * n.at(1, 1) + dl2 *n.at(1, 2) * n.at(1, 2);
            dL1.at(1, 2) = dl1 * n.at(2, 1) * n.at(2, 1) + dl2 *n.at(2, 2) * n.at(2, 2);
            dL1.at(1, 3) = dl1 * n.at(1, 1) * n.at(2, 1) + dl2 *n.at(1, 2) * n.at(2, 2);

            DB.beProductOf(D, B);
            LDB.beProductOf(dL1, DB);
            BktM22.beProductOf(BkT, M22);
            //BktM22Gk.beProductOf(BktM22,Gk);
            result1.beProductOf(BktM22Gk, LDB);

            // M12 + M21  = n1 \otimes n2 + n2 \otimes n1
            M12.at(1, 1) = n.at(1, 1) * n.at(1, 2) + n.at(1, 2) * n.at(1, 1);
            M12.at(1, 2) = n.at(1, 1) * n.at(2, 2) + n.at(1, 2) * n.at(2, 1);
            M12.at(2, 1) = n.at(2, 1) * n.at(1, 2) + n.at(2, 2) * n.at(1, 1);
            M12.at(2, 2) = n.at(2, 1) * n.at(2, 2) + n.at(2, 2) * n.at(2, 1);
            //dL2
            dL2.at(1, 1) = dl3 * ( n.at(1, 1) * n.at(1, 2) + n.at(1, 1) * n.at(1, 2) );
            dL2.at(1, 2) = dl3 * ( n.at(2, 1) * n.at(2, 2) + n.at(2, 1) * n.at(2, 2) );
            dL2.at(1, 3) = dl3 * ( n.at(1, 2) * n.at(2, 1) + n.at(1, 1) * n.at(2, 2) );

            LDB.beProductOf(dL2, DB);
            BktM12.beProductOf(BkT, M12);
            //BktM12Gk.beProductOf(BktM12,Gk);
            result2.beProductOf(BktM12Gk, LDB);

            result = result1;
            result.add(result2);

            result.times(dV);
            answer.add(result);
        }
    }
}


void
PhaseFieldElement :: computeStiffnessMatrix_kk(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    StructuralElement *elem = this->giveStructuralElement();
    double dV;
    double l;
    FloatMatrix lStiff;
    FloatMatrix Bk, Bt, BtB, N, Nt, NtN;
    StructuralCrossSection *cs = this->giveNLStructuralElement()->giveStructuralCrossSection();

    IntegrationRule *iRule = elem->giveIntegrationRule(0);

    answer.resize(nlSize, nlSize);

    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(j);

        GradDpMaterialExtensionInterface *dpmat = dynamic_cast< GradDpMaterialExtensionInterface * >( 
                cs->giveMaterialInterface( GradDpMaterialExtensionInterfaceType, gp ) );
        if ( !dpmat ) {
            OOFEM_ERROR("PhaseFieldElement :: computeStiffnessMatrix_uk - Material doesn't implement the required DpGrad interface!");
        }

        this->computeNkappaMatrixAt(gp, N);
        Nt.beTranspositionOf(N);
        this->computeBkappaMatrixAt(gp, Bk);
        Bt.beTranspositionOf(Bk);
        dV = elem->computeVolumeAround(gp);
        dpmat->givePDGradMatrix_kk(lStiff, rMode, gp, tStep);
        NtN.beProductOf(Nt, N);
        NtN.times(dV);
        answer.add(NtN);
        if ( averType == 0 || averType == 1 ) {
            l = lStiff.at(1, 1);
            BtB.beProductOf(Bt, Bk);
            BtB.times(l * l * dV);
            answer.add(BtB);
        } else if ( averType == 2 ) {
            FloatMatrix BtL, BtLB;
            BtL.beProductOf(Bt, lStiff);
            BtLB.beProductOf(BtL, Bk);
            BtLB.times(dV);
            answer.add(BtLB);
        }
    }
}


void
PhaseFieldElement :: computeStiffnessMatrix_uk(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();
    double dV;
    StructuralCrossSection *cs = this->giveNLStructuralElement()->giveStructuralCrossSection();
  
    IntegrationRule *iRule = elem->giveIntegrationRule(0);
    FloatMatrix B, Bt, Nk, BSN, BS, gPSigma;

    answer.resize(locSize, nlSize);

    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(j);

        GradDpMaterialExtensionInterface *dpmat = dynamic_cast< GradDpMaterialExtensionInterface * >( 
                cs->giveMaterialInterface( GradDpMaterialExtensionInterfaceType, gp ) );
        if ( !dpmat ) {
            OOFEM_ERROR("PhaseFieldElement :: computeStiffnessMatrix_uk - Material doesn't implement the required DpGrad interface!");
        }
        dpmat->givePDGradMatrix_uk(gPSigma, rMode, gp, tStep);
        this->computeNkappaMatrixAt(gp, Nk);
        elem->computeBmatrixAt(gp, B);
        if ( nlGeo == 1 ) {
            if ( elem->domain->giveEngngModel()->giveFormulation() == AL ) {
                elem->computeBmatrixAt(gp, B);
            } else {
                elem->computeBHmatrixAt(gp, B);
            }
        }

        Bt.beTranspositionOf(B);
        dV = elem->computeVolumeAround(gp);
        BS.beProductOf(Bt, gPSigma);
        BSN.beProductOf(BS, Nk);
        BSN.times(-dV);
        answer.add(BSN);
    }
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