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

#include "rankinematgrad.h"
#include "stressvector.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "error.h"
#include "classfactory.h"

namespace oofem {
/////////////////////////////////////////////////
// gradient regularization of Rankine plasticity
// coupled with isotropic damage
/////////////////////////////////////////////////

REGISTER_Material(RankineMatGrad);

// constructor
RankineMatGrad :: RankineMatGrad(int n, Domain *d) : RankineMat(n, d), GradDpMaterialExtensionInterface(d)
{
    L = 0.;
    negligible_damage = 0.;
}

int
RankineMatGrad :: hasMaterialModeCapability(MaterialMode mode)
{
    return mode == _PlaneStress;
}

void
RankineMatGrad :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}

void
RankineMatGrad :: givePDGradMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _PlaneStress:
        givePlaneStressStiffMtrx(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}

void
RankineMatGrad :: givePDGradMatrix_uk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _PlaneStress:
        givePlaneStressGprime(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}

void
RankineMatGrad :: givePDGradMatrix_ku(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _PlaneStress:
        givePlaneStressKappaMatrix(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}



void
RankineMatGrad :: givePDGradMatrix_kk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _PlaneStress:
        giveInternalLength(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}

void
RankineMatGrad :: givePDGradMatrix_LD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    OOFEM_ERROR("mMode = %d not supported\n", mMode);
}

void
RankineMatGrad :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );
    double tempDamage = status->giveTempDamage();
    double damage = status->giveDamage();
    double gprime;
    // Note:
    // The approximate solution of Helmholtz equation can lead
    // to very small but nonzero nonlocal kappa at some points that
    // are actually elastic. If such small values are positive,
    // they lead to a very small but nonzero damage. If this is
    // interpreted as "loading", the tangent terms are activated,
    // but damage will not actually grow at such points and the
    // convergence rate is slowed down. It is better to consider
    // such points as elastic.
    if ( tempDamage - damage <= negligible_damage ) {
        gprime = 0.;
    } else {
        double nonlocalCumulatedStrain = status->giveNonlocalCumulatedStrain();
        double tempLocalCumulatedStrain = status->giveTempCumulativePlasticStrain();
        double overNonlocalCumulatedStrain = mParam * nonlocalCumulatedStrain + ( 1. - mParam ) * tempLocalCumulatedStrain;
        gprime = computeDamageParamPrime(overNonlocalCumulatedStrain);
        gprime *= ( 1. - mParam );
    }

    evaluatePlaneStressStiffMtrx(answer, mode, gp, tStep, gprime);
}

// derivative of kappa (result of stress return) wrt final strain
void
RankineMatGrad :: givePlaneStressKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );
    answer.resize(1, 3);
    answer.zero();
    if ( mode != TangentStiffness ) {
        return;
    }

    double tempKappa = status->giveTempCumulativePlasticStrain();
    double dKappa = tempKappa - status->giveCumulativePlasticStrain();
    if ( dKappa <= 0. ) {
        return;
    }

    FloatArray eta(3);
    double dkap1 = status->giveDKappa(1);
    double H = evalPlasticModulus(tempKappa);

    // evaluate in principal coordinates

    if ( dkap1 == 0. ) {
        // regular case
        double Estar = E / ( 1. - nu * nu );
        double aux = Estar / ( H + Estar );
        eta.at(1) = aux;
        eta.at(2) = nu * aux;
        eta.at(3) = 0.;
    } else {
        // vertex case
        double dkap2 = status->giveDKappa(2);
        double denom = E * dkap1 + H * ( 1. - nu ) * ( dkap1 + dkap2 );
        eta.at(1) = E * dkap1 / denom;
        eta.at(2) = E * dkap2 / denom;
        eta.at(3) = 0.;
    }

    // transform to global coordinates

    FloatArray sigPrinc(2);
    FloatMatrix nPrinc(2, 2);
    StressVector effStress(status->giveTempEffectiveStress(), _PlaneStress);
    effStress.computePrincipalValDir(sigPrinc, nPrinc);

    FloatMatrix T(3, 3);
    givePlaneStressVectorTranformationMtrx(T, nPrinc, true);
    FloatArray etaglob(3);
    etaglob.beProductOf(T, eta);

    answer.at(1, 1) = etaglob.at(1);
    answer.at(1, 2) = etaglob.at(2);
    answer.at(1, 3) = etaglob.at(3);
}

// minus derivative of total stress wrt nonlocal kappa
void
RankineMatGrad :: givePlaneStressGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(3, 1);
    answer.zero();
    if ( mode != TangentStiffness ) {
        return;
    }

    RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );
    double damage = status->giveDamage();
    double tempDamage = status->giveTempDamage();
    if ( tempDamage - damage <= negligible_damage ) {
        return;
    }

    double nonlocalCumulatedStrain = status->giveNonlocalCumulatedStrain();
    double tempCumulatedStrain = status->giveTempCumulativePlasticStrain();
    double overNonlocalCumulatedStrain = mParam * nonlocalCumulatedStrain + ( 1. - mParam ) * tempCumulatedStrain;
    const FloatArray &tempEffStress = status->giveTempEffectiveStress();
    answer.at(1, 1) = tempEffStress.at(1);
    answer.at(2, 1) = tempEffStress.at(2);
    answer.at(3, 1) = tempEffStress.at(3);
    double gPrime = computeDamageParamPrime(overNonlocalCumulatedStrain);
    answer.times(gPrime * mParam);
}

void
RankineMatGrad :: giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(1, 1);
    answer.at(1, 1) = L;
}

void
RankineMatGrad :: giveRealStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep)
{
    RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    double tempDamage;
    RankineMat :: performPlasticityReturn(gp, totalStrain);

    tempDamage = computeDamage(gp, tStep);
    const FloatArray &tempEffStress = status->giveTempEffectiveStress();
    answer1.beScaled(1.0 - tempDamage, tempEffStress);
    answer2 = status->giveTempCumulativePlasticStrain();

    status->setNonlocalCumulatedStrain(nonlocalCumulatedStrain);
    status->letTempStrainVectorBe(totalStrain);
    status->setTempDamage(tempDamage);
    status->letTempEffectiveStressBe(tempEffStress);
    status->letTempStressVectorBe(answer1);
#ifdef keep_track_of_dissipated_energy
    double gf = sig0 * sig0 / E; // only estimated, but OK for this purpose
    status->computeWork_PlaneStress(gp, gf);
#endif
    double knl = giveNonlocalCumPlasticStrain(gp);
    double khat = mParam * knl + ( 1. - mParam ) * answer2;
    status->setKappa_nl(knl);
    status->setKappa_hat(khat);
}


void
RankineMatGrad :: computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep)
{
    RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );
    double localCumPlastStrain = status->giveTempCumulativePlasticStrain();
    double nlCumPlastStrain = status->giveNonlocalCumulatedStrain();
    kappa = mParam * nlCumPlastStrain + ( 1. - mParam ) * localCumPlastStrain;
}

double
RankineMatGrad :: giveNonlocalCumPlasticStrain(GaussPoint *gp)
{
    RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );
    return status->giveNonlocalCumulatedStrain();
}

IRResultType
RankineMatGrad :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, L, _IFT_RankineMatGrad_L);
    if ( L < 0.0 ) {
        L = 0.0;
    }

    mParam = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, _IFT_RankineMatGrad_m);

    negligible_damage = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, negligible_damage, _IFT_RankineMatGrad_negligibleDamage);

    return RankineMat :: initializeFrom(ir);
}


int
RankineMatGrad :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_CumPlasticStrain_2 ) {
        answer.resize(1);
        answer.at(1) = giveNonlocalCumPlasticStrain(gp);
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        computeCumPlastStrain(answer.at(1), gp, tStep);
        return 1;
    } else {
        return RankineMat :: giveIPValue(answer, gp, type, tStep);
    }
}

//=============================================================================
// GRADIENT RANKINE MATERIAL STATUS
//=============================================================================


RankineMatGradStatus :: RankineMatGradStatus(int n, Domain *d, GaussPoint *g) :
    RankineMatStatus(n, d, g)
{
    nonlocalCumulatedStrain = 0;
}

void
RankineMatGradStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "status {");
    fprintf(file, "damage %g, kappa %g, kappa_nl %g, kappa_hat %g", damage, kappa, kappa_nl, kappa_hat);
#ifdef keep_track_of_dissipated_energy
    fprintf(file, ", dissW %g, freeE %g, stressW %g", this->dissWork, ( this->stressWork ) - ( this->dissWork ), this->stressWork);
#endif
    fprintf(file, " }\n");
}


void
RankineMatGradStatus :: initTempStatus()
{
    // RankineMatStatus::initTempStatus();
    StructuralMaterialStatus :: initTempStatus();

    if ( plasticStrain.giveSize() == 0 ) {
        if ( gp->giveMaterialMode() == _PlaneStress ) {
            plasticStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector(_PlaneStress) );
        } else if ( gp->giveMaterialMode() == _3dMat ) {
            plasticStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector(_3dMat) );
        }

        plasticStrain.zero();
    }

    tempDamage = damage;
    tempPlasticStrain = plasticStrain;
    tempKappa = kappa;

    kappa_nl = 0.; // only containers
    kappa_hat = 0.;
}


void
RankineMatGradStatus :: updateYourself(TimeStep *tStep)
{
    RankineMatStatus :: updateYourself(tStep);
}
} // end namespace oofem
