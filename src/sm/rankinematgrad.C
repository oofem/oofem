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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "rankinematgrad.h"
#include "stressvector.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "error.h"

namespace oofem {
/////////////////////////////////////////////////
// gradient regularization of Rankine plasticity
// coupled with isotropic damage
/////////////////////////////////////////////////

// constructor
RankineMatGrad :: RankineMatGrad(int n, Domain *d) : RankineMat(n, d)
{
    R = 0.;
    negligible_damage = 0.;
}

int
RankineMatGrad :: hasMaterialModeCapability(MaterialMode mode)
{
    return ( mode == _PlaneStressGrad );
}

void
RankineMatGrad :: giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime)
//
// Returns characteristic material matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _PlaneStressGrad:
        if ( form == PDGrad_uu ) {
            givePlaneStressStiffMtrx(answer, form, rMode, gp, atTime);
            break;
        } else if ( form == PDGrad_ku ) {
            givePlaneStressKappaMatrix(answer, form, rMode, gp, atTime);
            break;
        } else if ( form == PDGrad_uk ) {
            givePlaneStressGprime(answer, form, rMode, gp, atTime);
            break;
        } else if ( form == PDGrad_kk ) {
            giveInternalLength(answer, form, rMode, gp, atTime);
            break;
        }

    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
        return;
    }
}

void
RankineMatGrad :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    RankineMatStatus *status = ( RankineMatStatus * ) this->giveStatus(gp);
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
        double tempKappa = status->giveTempCumulativePlasticStrain();
        double nlKappa = status->giveTempStrainVector().at(4);
        double kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
        gprime = computeDamageParamPrime(kappa);
        gprime *= ( 1. - mParam );
    }

    evaluatePlaneStressStiffMtrx(answer, form, mode, gp, atTime, gprime);
}

// derivative of kappa (result of stress return) wrt final strain
void
RankineMatGrad :: givePlaneStressKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    RankineMatGradStatus *status = ( RankineMatGradStatus * ) this->giveStatus(gp);
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
    StressVector effStress(_PlaneStress);
    status->giveTempEffectiveStress(effStress);
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
RankineMatGrad :: givePlaneStressGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    answer.resize(3, 1);
    answer.zero();
    if ( mode != TangentStiffness ) {
        return;
    }

    RankineMatGradStatus *status = ( RankineMatGradStatus * ) this->giveStatus(gp);
    double damage = status->giveDamage();
    double tempDamage = status->giveTempDamage();
    if ( tempDamage - damage <= negligible_damage ) {
        return;
    }

    double tempKappa = status->giveTempCumulativePlasticStrain();
    double nlKappa =  status->giveTempStrainVector().at(4);
    double kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
    FloatArray tempEffStress;
    status->giveTempEffectiveStress(tempEffStress);
    answer.at(1, 1) = tempEffStress.at(1);
    answer.at(2, 1) = tempEffStress.at(2);
    answer.at(3, 1) = tempEffStress.at(3);
    double gPrime = computeDamageParamPrime(kappa);
    answer.times(gPrime * mParam);
}

void
RankineMatGrad :: giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    answer.resize(1, 1);
    answer.at(1, 1) = R;
}

void
RankineMatGrad :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                       const FloatArray &totalStrain, TimeStep *atTime)
{
    RankineMatGradStatus *status = ( RankineMatGradStatus * ) this->giveStatus(gp);
    this->initGpForNewStep(gp);
    this->initTempStatus(gp);
    MaterialMode mode = gp->giveMaterialMode();
    MaterialMode plReturnMode;
    if ( mode == _PlaneStressGrad ) {
        plReturnMode = _PlaneStress;
    } else {
        _error("Unknown mode in RankineMatGrad :: giveRealStressVector\n");
    }

    double tempDam;
    FloatArray tempEffStress, totalStress, locTotalStrain;

    int size = totalStrain.giveSize();
    locTotalStrain = totalStrain;
    locTotalStrain.resize(size - 1);
    RankineMat ::  performPlasticityReturn(gp, locTotalStrain, plReturnMode);
    status->letTempStrainVectorBe(totalStrain);
    double localCumPlastStrain = status->giveTempCumulativePlasticStrain();
    tempDam = computeDamage(gp, atTime);
    status->giveTempEffectiveStress(tempEffStress);
    answer.beScaled( 1.0 - tempDam, tempEffStress);
    size = tempEffStress.giveSize();
    answer.resize(size + 1);
    answer.at(size + 1) = localCumPlastStrain;

    status->setTempDamage(tempDam);
    status->letTempEffectiveStressBe(tempEffStress);
    status->letTempStressVectorBe(answer);
#ifdef keep_track_of_dissipated_energy
    double gf = sig0 * sig0 / E; // only estimated, but OK for this purpose
    status->computeWork(gp, mode, gf);
#endif
    double knl = giveNonlocalCumPlasticStrain(gp);
    double khat = mParam * knl + ( 1. - mParam ) * localCumPlastStrain;
    status->setKappa_nl(knl);
    status->setKappa_hat(khat);
}


void
RankineMatGrad :: computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *atTime)
{
    double nlCumPlastStrain;
    RankineMatGradStatus *status = ( RankineMatGradStatus * ) this->giveStatus(gp);
    double localCumPlastStrain = status->giveTempCumulativePlasticStrain();
    FloatArray strain;
    strain = status->giveTempStrainVector();
    int size = strain.giveSize();
    nlCumPlastStrain = strain.at(size);
    kappa = mParam * nlCumPlastStrain + ( 1. - mParam ) * localCumPlastStrain;
}

double
RankineMatGrad :: giveNonlocalCumPlasticStrain(GaussPoint *gp)
{
    RankineMatGradStatus *status = ( RankineMatGradStatus * ) this->giveStatus(gp);
    FloatArray strain = status->giveTempStrainVector();
    int size = strain.giveSize();
    double answer = strain.at(size);
    return answer;
}

IRResultType
RankineMatGrad :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    RankineMat :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, R, IFT_RankineMatGrad_r, "r"); // Macro
    if ( R < 0.0 ) {
        R = 0.0;
    }

    mParam = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, IFT_RankineMatGrad_m, "m"); // Macro

    negligible_damage = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, negligible_damage, IFT_RankineMatGrad_m, "negligible_damage");

    return IRRT_OK;
}


int
RankineMatGrad :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    if ( type == IST_CumPlasticStrain_2 ) {
        answer.resize(1);
        answer.at(1) = giveNonlocalCumPlasticStrain(aGaussPoint);
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        computeCumPlastStrain(answer.at(1), aGaussPoint, atTime);
        return 1;
    } else {
        return RankineMat :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}

InternalStateValueType
RankineMatGrad :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_CumPlasticStrain_2 || type == IST_MaxEquivalentStrainLevel ) {
        return ISVT_SCALAR;
    } else {
        return RankineMat :: giveIPValueType(type);
    }
}

int
RankineMatGrad :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( type == IST_CumPlasticStrain_2 || type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return RankineMat :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}

int
RankineMatGrad :: giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    if ( type == IST_CumPlasticStrain_2 || type == IST_MaxEquivalentStrainLevel ) {
        return 1;
    } else {
        return RankineMat :: giveIPValueSize(type, gp);
    }
}

//=============================================================================
// GRADIENT RANKINE MATERIAL STATUS
//=============================================================================


RankineMatGradStatus :: RankineMatGradStatus(int n, Domain *d, GaussPoint *g) :
    RankineMatStatus(n, d, g)
{}

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
        if ( gp->giveMaterialMode() == _PlaneStressGrad ) {
            plasticStrain.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->giveSizeOfReducedStressStrainVector(_PlaneStress) );
        } else if ( gp->giveMaterialMode() == _3dMatGrad ) {
            plasticStrain.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->giveSizeOfReducedStressStrainVector(_3dMat) );
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
RankineMatGradStatus :: updateYourself(TimeStep *atTime)
{
    RankineMatStatus :: updateYourself(atTime);
}
} // end namespace oofem
