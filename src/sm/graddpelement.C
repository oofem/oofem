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


#include "graddpelement.h"
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
GradDpElement :: GradDpElement()
// Constructor.
{
    nlGeo = 0;
    averType = 0;
}

void
GradDpElement :: setDisplacementLocationArray(IntArray &answer, int nPrimNodes, int nPrimVars, int nSecNodes, int nSecVars)
{
    answer.resize(locSize);

    for ( int i = 1; i <= totalSize; i++ ) {
        if ( i < nSecNodes * nPrimVars + 1 ) {
            answer.at(i) = i + ( int ) ( ( ( i - 1 ) / nPrimVars ) ) * nSecVars;
        } else if ( i > nSecNodes * ( nPrimVars + nSecVars ) ) {
            answer.at(i - nSecVars * nSecNodes) = i;
        }
    }
}

void
GradDpElement :: setNonlocalLocationArray(IntArray &answer, int nPrimNodes, int nPrimVars, int nSecNodes, int nSecVars)
{
    answer.resize(nlSize);
    for ( int i = 1; i <= nlSize; i++ ) {
        answer.at(i) = i * nPrimVars + i;
    }
}


void
GradDpElement :: computeDisplacementDegreesOfFreedom(FloatArray &answer, TimeStep *tStep)
{
    StructuralElement *elem = this->giveStructuralElement();
    FloatArray u;
    answer.resize(locSize);

    elem->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
    u.resizeWithValues(totalSize);
    for ( int i = 1; i <= locSize; i++ ) {
        answer.at(i) = u.at( locU.at(i) );
    }
}

void GradDpElement :: computeNonlocalDegreesOfFreedom(FloatArray &answer, TimeStep *tStep)
{
    StructuralElement *elem = this->giveStructuralElement();
    FloatArray u;
    answer.resize(nlSize);

    elem->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
    u.resizeWithValues(totalSize);
    for  ( int i = 1; i <= nlSize; i++ ) {
        answer.at(i) = u.at( locK.at(i) );
    }
}

void
GradDpElement :: computeStressVectorAndLocalCumulatedStrain(FloatArray &answer, double localCumulatedStrain, GaussPoint *gp, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();

    double nlCumulatedStrain;

    nlGeo = elem->giveGeometryMode();
    StructuralCrossSection *cs = this->giveNLStructuralElement()->giveStructuralCrossSection();
    GradDpMaterialExtensionInterface *dpmat = static_cast< GradDpMaterialExtensionInterface * >(
        cs->giveMaterialInterface(GradDpMaterialExtensionInterfaceType, gp) );

    if ( !dpmat ) {
        OOFEM_ERROR("GradDpElement :: computeStiffnessMatrix_uu - Material doesn't implement the required DpGrad interface!");
    }

    this->computeNonlocalCumulatedStrain(nlCumulatedStrain, gp, tStep);
    if ( nlGeo == 0 ) {
        FloatArray Epsilon;
        this->computeLocalStrainVector(Epsilon, gp, tStep);
        if ( cs->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
            dpmat->giveRealStressVectorGrad(answer, localCumulatedStrain, gp, Epsilon, nlCumulatedStrain, tStep);
            return;
        } else {
            OOFEM_ERROR("computeStressVectorAndLocalCumulatedStrain : unsupported mode");
        }
    } else if ( nlGeo == 1 ) {
        if ( elem->giveDomain()->giveEngngModel()->giveFormulation() == TL ) {
            FloatArray vF;
            this->computeDeformationGradientVector(vF, gp, tStep);
            if ( cs->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
                dpmat->giveFirstPKStressVectorGrad(answer, localCumulatedStrain, gp, vF, nlCumulatedStrain, tStep);
                return;
            } else {
                OOFEM_ERROR("computeStressVectorAndLocalCumulatedStrain : unsupported mode");
            }
        } else {
            FloatArray vF;
            this->computeDeformationGradientVector(vF, gp, tStep);
            if ( cs->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
                dpmat->giveCauchyStressVectorGrad(answer, localCumulatedStrain, gp, vF, nlCumulatedStrain, tStep);
                return;
            } else {
                OOFEM_ERROR("computeStressVectorAndLocalCumulatedStrain : unsupported mode");
            }
        }
    }
}


void
GradDpElement :: computeLocalStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray u;
    FloatMatrix b;
    NLStructuralElement *elem = this->giveNLStructuralElement();
    nlGeo = elem->giveGeometryMode();

    this->computeDisplacementDegreesOfFreedom(u, tStep);
    elem->computeBmatrixAt(gp, b);
    answer.beProductOf(b, u);
}

void
GradDpElement :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    FloatMatrix B;
    NLStructuralElement *elem = this->giveNLStructuralElement();

    this->computeDisplacementDegreesOfFreedom(u, tStep);
    // Displacement gradient H = du/dX
    elem->computeBHmatrixAt(gp, B);
    answer.beProductOf(B, u);

    // Deformation gradient F = H + I
    MaterialMode matMode = gp->giveMaterialMode();
    if ( matMode == _3dMat || matMode == _PlaneStrain ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
        answer.at(3) += 1.0;
    } else if ( matMode == _PlaneStress ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
    } else if ( matMode == _1dMat ) {
        answer.at(1) += 1.0;
    } else {
        OOFEM_ERROR2( "computeDeformationGradientVector : MaterialMode is not supported yet (%s)", __MaterialModeToString(matMode) );
    }
}

void
GradDpElement :: computeNonlocalCumulatedStrain(double &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix Nk;
    FloatArray u;
    FloatArray aux;

    this->computeNkappaMatrixAt(gp, Nk);
    this->computeNonlocalDegreesOfFreedom(u, tStep);
    aux.beProductOf(Nk, u);
    answer = aux.at(1);
}


void
GradDpElement :: computeNonlocalGradient(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix Bk;
    FloatArray u;
    FloatArray aux;

    this->computeBkappaMatrixAt(gp, Bk);
    this->computeNonlocalDegreesOfFreedom(u, tStep);
    answer.beProductOf(Bk, u);
}


void
GradDpElement :: giveNonlocalInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    double dV, localCumulatedStrain = 0.;
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
GradDpElement :: giveLocalInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();
    IntegrationRule *iRule = elem->giveIntegrationRule(0);
    nlGeo = elem->giveGeometryMode();
    FloatArray BS, vStress;
    FloatMatrix B;
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);

        if ( nlGeo == 0 || elem->domain->giveEngngModel()->giveFormulation() == AL ) {
            elem->computeBmatrixAt(gp, B);
        } else if ( nlGeo == 1 ) {
            elem->computeBHmatrixAt(gp, B);
        }
        vStress = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveTempStressVector();

        if ( vStress.giveSize() == 0 ) {         /// @todo is this check really necessary?
            break;
        }

        // Compute nodal internal forces at nodes as f = B^T*Stress dV
        double dV  = elem->computeVolumeAround(gp);
        BS.beTProductOf(B, vStress);
        answer.add(dV, BS);
    }
}

void
GradDpElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
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
GradDpElement :: computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    //set displacement and nonlocal location array
    this->setDisplacementLocationArray(locU, nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    this->setNonlocalLocationArray(locK, nPrimNodes, nPrimVars, nSecNodes, nSecVars);


    FloatArray localForces(locSize);
    FloatArray nlForces(nlSize);
    answer.resize(totalSize);

    this->computeLocForceLoadVector(localForces, tStep, mode);

    answer.assemble(localForces, locU);
    answer.assemble(nlForces, locK);
}


/************************************************************************/
void
GradDpElement :: computeLocForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further sobstract part corresponding to non-nodeal loading.
{
    FloatMatrix T;
    NLStructuralElement *elem = this->giveNLStructuralElement();
    elem->computeLocalForceLoadVector(answer, tStep, mode);

    // transform result from global cs to nodal cs. if necessary
    if ( answer.isNotEmpty() ) {
        if ( elem->computeGtoLRotationMatrix(T) ) {
            // first back to global cs from element local
            answer.rotatedWith(T, 't');
        }
    } else {
        answer.resize(locSize);
        answer.zero();
    }
}


void
GradDpElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
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
GradDpElement :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
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
            cs->giveMaterialInterface(GradDpMaterialExtensionInterfaceType, gp) );
        if ( !dpmat ) {
            OOFEM_ERROR("GradDpElement :: computeStiffnessMatrix_uk - Material doesn't implement the required DpGrad interface!");
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
GradDpElement :: computeStiffnessMatrix_ku(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
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
            cs->giveMaterialInterface(GradDpMaterialExtensionInterfaceType, gp) );
        if ( !dpmat ) {
            OOFEM_ERROR("GradDpElement :: computeStiffnessMatrix_uk - Material doesn't implement the required DpGrad interface!");
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
GradDpElement :: computeStiffnessMatrix_kk(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
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
            cs->giveMaterialInterface(GradDpMaterialExtensionInterfaceType, gp) );
        if ( !dpmat ) {
            OOFEM_ERROR("GradDpElement :: computeStiffnessMatrix_uk - Material doesn't implement the required DpGrad interface!");
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
GradDpElement :: computeStiffnessMatrix_uk(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
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
            cs->giveMaterialInterface(GradDpMaterialExtensionInterfaceType, gp) );
        if ( !dpmat ) {
            OOFEM_ERROR("GradDpElement :: computeStiffnessMatrix_uk - Material doesn't implement the required DpGrad interface!");
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
GradDpElement :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                // Required by IR_GIVE_FIELD macro
    //nlGeo = 0;

    return IRRT_OK;
}
} // end namespace oofem
