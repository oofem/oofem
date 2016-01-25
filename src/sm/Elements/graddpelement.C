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


#include "../sm/Elements/graddpelement.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/Elements/nlstructuralelement.h"
#include "../sm/Materials/graddpmaterialextensioninterface.h"
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
#include "nonlocalbarrier.h"
#include "engngm.h"

#include <cstdio>

namespace oofem {
GradDpElement :: GradDpElement()
{ }

void
GradDpElement :: setDisplacementLocationArray()
{
    locU.resize(locSize);

    for ( int i = 1; i <= totalSize; i++ ) {
        if ( i < nSecNodes * nPrimVars + 1 ) {
            locU.at(i) = i + ( int ) ( ( ( i - 1 ) / nPrimVars ) ) * nSecVars;
        } else if ( i > nSecNodes * ( nPrimVars + nSecVars ) ) {
            locU.at(i - nSecVars * nSecNodes) = i;
        }
    }
}

void
GradDpElement :: setNonlocalLocationArray()
{
    locK.resize(nlSize);
    for ( int i = 1; i <= nlSize; i++ ) {
        locK.at(i) = i * nPrimVars + i;
    }
}


void
GradDpElement :: computeDisplacementDegreesOfFreedom(FloatArray &answer, TimeStep *tStep)
{
    this->giveStructuralElement()->computeVectorOf({D_u, D_v, D_w}, VM_Total, tStep, answer);
}

void GradDpElement :: computeNonlocalDegreesOfFreedom(FloatArray &answer, TimeStep *tStep)
{
    this->giveStructuralElement()->computeVectorOf({G_0}, VM_Total, tStep, answer);
}

void
GradDpElement :: computeStressVectorAndLocalCumulatedStrain(FloatArray &answer, double localCumulatedStrain, GaussPoint *gp, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();

    double nlCumulatedStrain;

    int nlGeo = elem->giveGeometryMode();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    GradDpMaterialExtensionInterface *dpmat = static_cast< GradDpMaterialExtensionInterface * >(
        cs->giveMaterialInterface(GradDpMaterialExtensionInterfaceType, gp) );

    if ( !dpmat ) {
        OOFEM_ERROR("Material doesn't implement the required DpGrad interface!");
    }

    this->computeNonlocalCumulatedStrain(nlCumulatedStrain, gp, tStep);
    if ( nlGeo == 0 ) {
        FloatArray Epsilon;
        this->computeLocalStrainVector(Epsilon, gp, tStep);
        dpmat->giveRealStressVectorGrad(answer, localCumulatedStrain, gp, Epsilon, nlCumulatedStrain, tStep);
    } else if ( nlGeo == 1 ) {
        if ( elem->giveDomain()->giveEngngModel()->giveFormulation() == TL ) {
            FloatArray vF;
            this->computeDeformationGradientVector(vF, gp, tStep);
            dpmat->giveFirstPKStressVectorGrad(answer, localCumulatedStrain, gp, vF, nlCumulatedStrain, tStep);
        } else {
            FloatArray vF;
            this->computeDeformationGradientVector(vF, gp, tStep);
            dpmat->giveCauchyStressVectorGrad(answer, localCumulatedStrain, gp, vF, nlCumulatedStrain, tStep);
        }
    }
}


void
GradDpElement :: computeLocalStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray u;
    FloatMatrix b;
    NLStructuralElement *elem = this->giveNLStructuralElement();

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
        OOFEM_ERROR( "MaterialMode is not supported yet (%s)", __MaterialModeToString(matMode) );
    }
}

void
GradDpElement :: computeNonlocalCumulatedStrain(double &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray Nk, u;

    this->computeNkappaMatrixAt(gp, Nk);
    this->computeNonlocalDegreesOfFreedom(u, tStep);
    answer = Nk.dotProduct(u);
}


void
GradDpElement :: computeNonlocalGradient(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix Bk;
    FloatArray u;

    this->computeBkappaMatrixAt(gp, Bk);
    this->computeNonlocalDegreesOfFreedom(u, tStep);
    answer.beProductOf(Bk, u);
}


void
GradDpElement :: giveNonlocalInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    double localCumulatedStrain = 0.;
    NLStructuralElement *elem = this->giveNLStructuralElement();
    FloatMatrix stiffKappa;
    FloatArray Nk;
    FloatArray aux, dKappa, stress;

    int size = nSecVars * nSecNodes;

    //set displacement and nonlocal location array
    this->setDisplacementLocationArray();
    this->setNonlocalLocationArray();

    answer.resize(size);
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
        this->computeNkappaMatrixAt(gp, Nk);

        double dV = elem->computeVolumeAround(gp);
        this->computeStressVectorAndLocalCumulatedStrain(stress, localCumulatedStrain, gp, tStep);
        aux.add(-dV * localCumulatedStrain, Nk);
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
    int nlGeo = elem->giveGeometryMode();
    FloatArray BS, vStress;
    FloatMatrix B;
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {

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

#if 0
void
GradDpElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatArray answerU, answerK;

    double localCumulatedStrain = 0.;
    FloatMatrix stiffKappa, B;
    FloatArray Nk, aux, dKappa, stress;
    FloatMatrix lStiff;
    FloatMatrix Bk;
    FloatArray gKappa, L_gKappa;

    //set displacement and nonlocal location array
    this->setDisplacementLocationArray();
    this->setNonlocalLocationArray();

    this->computeNonlocalDegreesOfFreedom(dKappa, tStep);

    int nlGeo = elem->giveGeometryMode();
    for ( auto &gp: *elem->giveIntegrationRule(0) ) {
        GradDpMaterialExtensionInterface *dpmat = dynamic_cast< GradDpMaterialExtensionInterface * >(
            cs->giveMaterialInterface(GradDpMaterialExtensionInterfaceType, gp) );
        if ( !dpmat ) {
            OOFEM_ERROR("Material doesn't implement the required DpGrad interface!");
        }

        double dV = elem->computeVolumeAround(gp);

        if ( nlGeo == 0 || elem->domain->giveEngngModel()->giveFormulation() == AL ) {
            elem->computeBmatrixAt(gp, B);
        } else if ( nlGeo == 1 ) {
            elem->computeBHmatrixAt(gp, B);
        }
        this->computeStressVectorAndLocalCumulatedStrain(stress, localCumulatedStrain, gp, tStep);

        answerU.plusProduct(B, stress, dV);

        // Gradient part:
        this->computeNkappaMatrixAt(gp, Nk);
        this->computeBkappaMatrixAt(gp, Bk);

        dpmat->givePDGradMatrix_kk(lStiff, TangentStiffness, gp, tStep);
        double kappa = Nk.dotProduct(dKappa);
        gKappa.beProductOf(Bk, dKappa);

        answerK.add(-dV * localCumulatedStrain, Nk);
        answerK.add(kappa * dV, Nk);

        if ( dpmat->giveAveragingType() == 0 || dpmat->giveAveragingType() == 1 ) {
            double l = lStiff.at(1, 1);
            answerK.plusProduct(Bk, gKappa, l * l * dV);
        } else if ( dpmat->giveAveragingType() == 2 ) {
            L_gKappa.beProductOf(lStiff, gKappa);
            answerK.plusProduct(Bk, L_gKappa, dV);
        }
    }

    //this->computeStiffnessMatrix_kk(stiffKappa, TangentStiffness, tStep);
    //answerK.beProductOf(stiffKappa, dKappa);
    answerK.add(aux);

    answer.resize(totalSize);
    answer.zero();
    answer.assemble(answerU, locU);
    answer.assemble(answerK, locK);
}
#else
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

    ///@todo This code relies on locU and locK already exists. Design seems highly unsafe.
    /// Instead locU and locK shoulld be statically determined (cf. with the Taylor-Hood elements) / Mikael
    this->giveNonlocalInternalForcesVector(answerK, tStep, useUpdatedGpRecord);
    this->giveLocalInternalForcesVector(answerU, tStep, useUpdatedGpRecord);


    answer.assemble(answerU, locU);
    answer.assemble(answerK, locK);
}
#endif

void
GradDpElement :: computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    //set displacement and nonlocal location array
    this->setDisplacementLocationArray();
    this->setNonlocalLocationArray();


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

#if 0
void
GradDpElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
    this->setDisplacementLocationArray();
    this->setNonlocalLocationArray();

    NLStructuralElement *elem = this->giveNLStructuralElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatMatrix B, D, DB;
    FloatMatrix DkuB, Dku;
    FloatArray Nk;
    FloatMatrix SNk, gPSigma;
    FloatMatrix lStiff;
    FloatMatrix Bk, LBk;
    FloatMatrix answer_uu, answer_ku, answer_uk, answer_kk;

    int nlGeo = elem->giveGeometryMode();
    bool matStiffSymmFlag = elem->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);

    for ( auto &gp : *elem->giveIntegrationRule(0) ) {
        GradDpMaterialExtensionInterface *dpmat = dynamic_cast< GradDpMaterialExtensionInterface * >(
            cs->giveMaterialInterface(GradDpMaterialExtensionInterfaceType, gp) );
        if ( !dpmat ) {
            OOFEM_ERROR("Material doesn't implement the required DpGrad interface!");
        }

        double dV = elem->computeVolumeAround(gp);

        if ( nlGeo == 0 ) {
            elem->computeBmatrixAt(gp, B);
        } else if ( nlGeo == 1 ) {
            if ( elem->domain->giveEngngModel()->giveFormulation() == AL ) {
                elem->computeBmatrixAt(gp, B);
            } else {
                elem->computeBHmatrixAt(gp, B);
            }
        }
        this->computeNkappaMatrixAt(gp, Nk);
        this->computeBkappaMatrixAt(gp, Bk);

        dpmat->givePDGradMatrix_uu(D, rMode, gp, tStep);
        dpmat->givePDGradMatrix_ku(Dku, rMode, gp, tStep);
        dpmat->givePDGradMatrix_uk(gPSigma, rMode, gp, tStep);
        dpmat->givePDGradMatrix_kk(lStiff, rMode, gp, tStep);

        /////////////////////////////////////////////////////////////////// uu:
        DB.beProductOf(D, B);
        if ( matStiffSymmFlag ) {
            answer_uu.plusProductSymmUpper(B, DB, dV);
        } else {
            answer_uu.plusProductUnsym(B, DB, dV);
        }

        //////////////////////////////////////////////////////////////////////// ku:
        DkuB.beProductOf(Dku, B);
        answer_ku.plusProductUnsym(Nk, DkuB, -dV);

        if ( dpmat->giveAveragingType() == 2 ) {
            double dl1, dl2, dl3;
            FloatMatrix LDB;
            FloatMatrix GkLDB, MGkLDB;
            FloatMatrix M22, M12;
            FloatMatrix dL1(1, 3), dL2(1, 3), dLdS;
            FloatArray Gk, n1, n2;


            dpmat->givePDGradMatrix_LD(dLdS, rMode, gp, tStep);
            this->computeNonlocalGradient(Gk, gp, tStep);

            dl1 = dLdS.at(3, 3);
            dl2 = dLdS.at(4, 4);
            dl3 = dLdS.at(5, 5);
            n1 = {dLdS.at(1, 1), dLdS.at(2, 1)};
            n2 = {dLdS.at(1, 2), dLdS.at(2, 2)};

            // first term Bk^T M22 G L1 D B
            // M22 = n2 \otimes n2
            M22.plusDyadUnsym(n2, n2, 1.);
            // dL1
            dL1.at(1, 1) = dl1 * n1.at(1) * n1.at(1) + dl2 * n2.at(1) * n2.at(1);
            dL1.at(1, 2) = dl1 * n1.at(2) * n1.at(2) + dl2 * n2.at(2) * n2.at(2);
            dL1.at(1, 3) = dl1 * n1.at(1) * n1.at(2) + dl2 * n2.at(1) * n2.at(2);

            LDB.beProductOf(dL1, DB);
            GkLDB.beProductOf(Gk, LDB);
            MGkLDB.beProductOf(M22, GkLDB);
            answer.plusProductUnsym(Bk, MGkLDB, dV);

            // M12 + M21  = n1 \otimes n2 + n2 \otimes n1
            M12.plusDyadUnsym(n1, n2, 1.);
            M12.plusDyadUnsym(n2, n1, 1.);
            //dL2
            dL2.at(1, 1) = dl3 * ( n1.at(1) * n2.at(1) + n1.at(1) * n2.at(1) );
            dL2.at(1, 2) = dl3 * ( n1.at(2) * n2.at(2) + n1.at(2) * n2.at(2) );
            dL2.at(1, 3) = dl3 * ( n1.at(2) * n2.at(1) + n1.at(1) * n2.at(2) );

            // Bk * ((M12 * L2 + M22 * L1) * DB)
            LDB.beProductOf(dL2, DB);
            GkLDB.beProductOf(Gk, LDB);
            MGkLDB.beProductOf(M12, GkLDB);
            answer.plusProductUnsym(Bk, MGkLDB, dV);
        }

        //////////////////////////////////////////////////////////////////////// uk:
        SNk.beProductOf(gPSigma, Nk);
        answer_uk.plusProductUnsym(B, SNk, -dV); // uk

        /////////////////////////////////////////////////////////////////////// kk:
        answer_kk.plusProductUnsym(Nk, Nk, dV);
        if ( dpmat->giveAveragingType() == 0 || dpmat->giveAveragingType() == 1 ) {
            double l = lStiff.at(1, 1);
            answer_kk.plusProductUnsym(Bk, Bk, l * l * dV);
        } else if ( dpmat->giveAveragingType() == 2 ) {
            LBk.beProductOf(lStiff, Bk);
            answer_kk.plusProductUnsym(Bk, LBk, dV);
        }
    }

    if ( elem->domain->giveEngngModel()->giveFormulation() == AL ) {
        FloatMatrix initialStressMatrix;
        elem->computeInitialStressMatrix(initialStressMatrix, tStep);
        answer_uu.add(initialStressMatrix);
    }

    if ( matStiffSymmFlag ) {
        answer_uu.symmetrized();
    }

    answer.resize(totalSize, totalSize);
    answer.zero();
    answer.assemble(answer_uu, locU);
    answer.assemble(answer_uk, locU, locK);
    answer.assemble(answer_ku, locK, locU);
    answer.assemble(answer_kk,locK);
}
#else

void
GradDpElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
    this->setDisplacementLocationArray();
    this->setNonlocalLocationArray();

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
#endif

void
GradDpElement :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatMatrix B, D, DB;

    int nlGeo = elem->giveGeometryMode();
    bool matStiffSymmFlag = elem->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answer.clear();
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {

        GradDpMaterialExtensionInterface *dpmat = dynamic_cast< GradDpMaterialExtensionInterface * >(
            cs->giveMaterialInterface(GradDpMaterialExtensionInterfaceType, gp) );
        if ( !dpmat ) {
            OOFEM_ERROR("Material doesn't implement the required DpGrad interface!");
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
    FloatArray Nk;
    FloatMatrix B, DkuB, Dku;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();

    answer.clear();

    int nlGeo = elem->giveGeometryMode();

    for ( auto &gp: *elem->giveIntegrationRule(0) ) {

        GradDpMaterialExtensionInterface *dpmat = dynamic_cast< GradDpMaterialExtensionInterface * >(
            cs->giveMaterialInterface(GradDpMaterialExtensionInterfaceType, gp) );
        if ( !dpmat ) {
            OOFEM_ERROR("Material doesn't implement the required DpGrad interface!");
        }

        elem->computeBmatrixAt(gp, B);
        if ( nlGeo == 1 ) {
            if ( elem->domain->giveEngngModel()->giveFormulation() == AL ) {
                elem->computeBmatrixAt(gp, B);
            } else {
                elem->computeBHmatrixAt(gp, B);
            }
        }

        dpmat->givePDGradMatrix_ku(Dku, rMode, gp, tStep);
        this->computeNkappaMatrixAt(gp, Nk);
        dV = elem->computeVolumeAround(gp);
        DkuB.beProductOf(Dku, B);
        answer.plusProductUnsym(Nk, DkuB, -dV);

        if ( dpmat->giveAveragingType() == 2 ) {
            double dl1, dl2, dl3;
            FloatArray Gk;
            FloatMatrix D, DB, LDB;
            FloatMatrix Bk, BktM22, BktM22Gk, BktM12, BktM12Gk, M22(2, 2), M12(2, 2);
            FloatMatrix dL1(1, 3), dL2(1, 3), result1, result2, dLdS, n(2, 2);

            this->computeBkappaMatrixAt(gp, Bk);
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
            // first term Bk^T M22 G L1 D B
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
            BktM22.beTProductOf(Bk, M22);
            ///@todo This can't possibly work if this is uncommented (!) / Mikael
            //BktM22Gk.beProductOf(BktM22,Gk);
            result1.beProductOf(BktM22Gk, LDB);
            answer.add(dV, result1);
            // This would be slightly shorter and faster;
            //GkLDB.beProductOf(Gk, LDB);
            //MGkLDB.beProductOf(M22, GkLDB);
            //answer.plusProductUnsym(Bk, MGkLDB, dV);

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
            BktM12.beTProductOf(Bk, M12);
            ///@todo This can't possibly work if this is uncommented (!) / Mikael
            //BktM12Gk.beProductOf(BktM12,Gk);
            result2.beProductOf(BktM12Gk, LDB);
            answer.add(dV, result2);
            // This would be slightly shorter and faster;
            //GkLDB.beProductOf(Gk, LDB);
            //MGkLDB.beProductOf(M12, GkLDB);
            //answer.plusProductUnsym(Bk, MGkLDB, dV);
        }
    }
}


void
GradDpElement :: computeStiffnessMatrix_kk(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    StructuralElement *elem = this->giveStructuralElement();
    double dV;
    FloatMatrix lStiff;
    FloatArray Nk;
    FloatMatrix Bk, LBk;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();

    answer.clear();

    for ( auto &gp: *elem->giveIntegrationRule(0) ) {

        GradDpMaterialExtensionInterface *dpmat = dynamic_cast< GradDpMaterialExtensionInterface * >(
            cs->giveMaterialInterface(GradDpMaterialExtensionInterfaceType, gp) );
        if ( !dpmat ) {
            OOFEM_ERROR("Material doesn't implement the required DpGrad interface!");
        }

        this->computeNkappaMatrixAt(gp, Nk);
        this->computeBkappaMatrixAt(gp, Bk);
        dV = elem->computeVolumeAround(gp);

        dpmat->givePDGradMatrix_kk(lStiff, rMode, gp, tStep);
        answer.plusProductUnsym(Nk, Nk, dV);
        if ( dpmat->giveAveragingType() == 0 || dpmat->giveAveragingType() == 1 ) {
            double l = lStiff.at(1, 1);
            answer.plusProductUnsym(Bk, Bk, l * l * dV);
        } else if ( dpmat->giveAveragingType() == 2 ) {
            LBk.beProductOf(lStiff, Bk);
            answer.plusProductUnsym(Bk, LBk, dV);
        }
    }
}


void
GradDpElement :: computeStiffnessMatrix_uk(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();
    double dV;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    int nlGeo = elem->giveGeometryMode();
    FloatArray Nk;
    FloatMatrix B, SNk, gPSigma;

    answer.clear();

    for ( auto &gp: *elem->giveIntegrationRule(0) ) {

        GradDpMaterialExtensionInterface *dpmat = dynamic_cast< GradDpMaterialExtensionInterface * >(
            cs->giveMaterialInterface(GradDpMaterialExtensionInterfaceType, gp) );
        if ( !dpmat ) {
            OOFEM_ERROR("Material doesn't implement the required DpGrad interface!");
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
        dV = elem->computeVolumeAround(gp);

        SNk.beProductOf(gPSigma, Nk);
        answer.plusProductUnsym(B, SNk, -dV);
    }
}


IRResultType
GradDpElement :: initializeFrom(InputRecord *ir)
{
    //IRResultType result;                // Required by IR_GIVE_FIELD macro
    //nlGeo = 0;

    return IRRT_OK;
}
} // end namespace oofem
