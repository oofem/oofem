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

RankineMatGrad :: RankineMatGrad(int n, Domain *d) : RankineMat(n, d), GradientDamageMaterialExtensionInterface(d)
{}

/////////////////////////////////////////////////////////////////////////////
void
RankineMatGrad :: initializeFrom(InputRecord &ir)
{
    RankineMat :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, L, _IFT_RankineMatGrad_L);
    if ( L < 0.0 ) {
        L = 0.0;
    }

    mParam = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, _IFT_RankineMatGrad_m);

    negligible_damage = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, negligible_damage, _IFT_RankineMatGrad_negligibleDamage);


    int formulationType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, formulationType, _IFT_RankineMatGrad_formulationType);
    if ( formulationType == 0 ) {
        this->gradientDamageFormulationType = GDFT_Standard;
    } else if ( formulationType == 2 ) {
        this->gradientDamageFormulationType =   GDFT_Eikonal;
    } else {
        throw ValueInputException(ir, _IFT_RankineMatGrad_formulationType, "Unknown gradient damage formulation");
    }
}
/////////////////////////////////////////////////////////////////////////////

bool
RankineMatGrad :: hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _PlaneStress;
}

void
RankineMatGrad :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const 
//
// Returns characteristic material matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}

void
RankineMatGrad :: giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _PlaneStress:
        if ( mode == ElasticStiffness ) {
            this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
        } else {
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

            answer = evaluatePlaneStressStiffMtrx(mode, gp, tStep, gprime);
        }
        break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}

void
RankineMatGrad :: giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _PlaneStress:
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
        answer.times(-1. * gPrime * mParam);
    }
    break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}

void
RankineMatGrad :: giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _PlaneStress:
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

        FloatMatrix T = givePlaneStressVectorTranformationMtrx(nPrinc, true);
        FloatArray etaglob(3);
        etaglob.beProductOf(T, eta);

        answer.at(1, 1) = etaglob.at(1);
        answer.at(1, 2) = etaglob.at(2);
        answer.at(1, 3) = etaglob.at(3);

        if ( gradientDamageFormulationType == GDFT_Standard ) {
            answer.times(1.);
        } else if ( gradientDamageFormulationType == GDFT_Eikonal ) {
            double iA = this->computeEikonalInternalLength_a(gp);
            if ( iA != 0 ) {
                answer.times(1. / iA);
            }
        } else {
            OOFEM_WARNING("Unknown internalLengthDependenceType");
        }
    }
    break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}

void
RankineMatGrad :: giveGradientDamageStiffnessMatrix_du_NB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    if ( gradientDamageFormulationType == GDFT_Standard ) {
        answer.clear();
    } else if ( gradientDamageFormulationType == GDFT_Eikonal ) {
        MaterialMode mMode = gp->giveMaterialMode();
        switch ( mMode ) {
        case _PlaneStress:
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

            FloatMatrix T = givePlaneStressVectorTranformationMtrx(nPrinc, true);
            FloatArray etaglob(3);
            etaglob.beProductOf(T, eta);

            answer.at(1, 1) = etaglob.at(1);
            answer.at(1, 2) = etaglob.at(2);
            answer.at(1, 3) = etaglob.at(3);

            double LocalCumulatedStrain = status->giveTempLocalDamageDrivingVariable();
            double NonLocalCumulatedStrain = status->giveTempNonlocalDamageDrivingVariable();
            answer.times(LocalCumulatedStrain - NonLocalCumulatedStrain);
            double iA = this->computeEikonalInternalLength_a(gp);
            if ( iA != 0 ) {
                answer.times( 1. / ( iA * iA ) );
            }
            double iAPrime = this->computeEikonalInternalLength_aPrime(gp);
            double gPrime = this->computeDamageParamPrime(tempKappa);
            answer.times(iAPrime * gPrime);
            answer.times(1. - mParam);
        }
        break;
        default:
            OOFEM_ERROR("mMode = %d not supported\n", mMode);
        }
    } else {
        OOFEM_WARNING("Unknown internalLengthDependenceType");
    }
}

void
RankineMatGrad :: giveGradientDamageStiffnessMatrix_du_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    if ( gradientDamageFormulationType == GDFT_Standard ) {
        answer.clear();
    } else if ( gradientDamageFormulationType == GDFT_Eikonal )  {
        MaterialMode mMode = gp->giveMaterialMode();
        switch ( mMode ) {
        case _PlaneStress:
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

            FloatMatrix T = givePlaneStressVectorTranformationMtrx(nPrinc, true);
            FloatArray etaglob(3);
            etaglob.beProductOf(T, eta);

            FloatArray GradP = status->giveTempNonlocalDamageDrivingVariableGrad();
            answer.beDyadicProductOf(GradP, etaglob);
            double iBPrime = this->computeEikonalInternalLength_bPrime(gp);
            double gPrime = this->computeDamageParamPrime(tempKappa);
            answer.times( iBPrime * gPrime * ( 1. - mParam ) );
        }
        break;
        default:
            OOFEM_ERROR("mMode = %d not supported\n", mMode);
        }
    } else {
        OOFEM_WARNING("Unknown internalLengthDependenceType");
    }
}

void
RankineMatGrad :: giveGradientDamageStiffnessMatrix_dd_NN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _PlaneStress:
        if ( gradientDamageFormulationType == GDFT_Standard ) {
            answer.clear();
        } else if ( gradientDamageFormulationType == GDFT_Eikonal )  {
            RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );
            answer.resize(1, 1);
            answer.zero();
            double iA = this->computeEikonalInternalLength_a(gp);

            if ( iA != 0 ) {
                answer.at(1, 1) += 1. / iA;
            }
            if ( mode == TangentStiffness ) {
                double tempKappa = status->giveTempCumulativePlasticStrain();
                if ( tempKappa > status->giveCumulativePlasticStrain() && iA != 0 ) {
                    double iAPrime = this->computeEikonalInternalLength_aPrime(gp);
                    double gPrime = this->computeDamageParamPrime(tempKappa);
                    double LocalCumulatedStrain = status->giveTempLocalDamageDrivingVariable();
                    double NonLocalCumulatedStrain = status->giveTempNonlocalDamageDrivingVariable();
                    answer.at(1, 1) += iAPrime / iA / iA * gPrime * mParam * ( LocalCumulatedStrain - NonLocalCumulatedStrain );
                }
            }
        } else {
            OOFEM_WARNING("Unknown internalLengthDependenceType");
        }
        break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}

void
RankineMatGrad :: giveGradientDamageStiffnessMatrix_dd_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _PlaneStress:
    {
        int n = this->giveDimension(gp);
        answer.resize(n, n);
        answer.beUnitMatrix();
        if ( gradientDamageFormulationType == GDFT_Standard ) {
            answer.times(internalLength * internalLength);
        } else if ( gradientDamageFormulationType == GDFT_Eikonal ) {
            double iB = this->computeEikonalInternalLength_b(gp);
            answer.times(iB);
        } else {
            OOFEM_WARNING("Unknown internalLengthDependenceType");
        }
        break;
    }
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}

void
RankineMatGrad :: giveGradientDamageStiffnessMatrix_dd_BN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    if ( gradientDamageFormulationType == GDFT_Standard ) {
        answer.clear();
    } else if ( gradientDamageFormulationType == GDFT_Eikonal )  {
        MaterialMode mMode = gp->giveMaterialMode();
        switch ( mMode ) {
        case _PlaneStress:
        {
            RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );
            FloatArray GradP = status->giveTempNonlocalDamageDrivingVariableGrad();
            answer = GradP;
            double iBPrime = this->computeEikonalInternalLength_bPrime(gp);
            double tempKappa = status->giveTempCumulativePlasticStrain();
            double gPrime = this->computeDamageParamPrime(tempKappa);
            answer.times(iBPrime * gPrime * mParam);
            break;
        }
        default:
            OOFEM_ERROR("mMode = %d not supported\n", mMode);
        }
    } else {
        OOFEM_WARNING("Unknown internalLengthDependenceType");
    }
}

FloatMatrixF<3,3>
RankineMatGrad :: givePlaneStressStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
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

    return evaluatePlaneStressStiffMtrx(mode, gp, tStep, gprime);
}

void
RankineMatGrad :: giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(1, 1);
    answer.at(1, 1) = L;
}

void
RankineMatGrad :: giveNonlocalInternalForces_N_factor(double &answer, double nlDamageDrivingVariable, GaussPoint *gp, TimeStep *tStep)
{
    // I modified this one to put pnl -p instead of just pnl
    RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );
    double LocalCumulatedStrain = status->giveTempLocalDamageDrivingVariable();
    double NonLocalCumulatedStrain = status->giveTempNonlocalDamageDrivingVariable();
    answer = NonLocalCumulatedStrain - LocalCumulatedStrain;
    if ( gradientDamageFormulationType == GDFT_Eikonal ) {
        double iA = this->computeEikonalInternalLength_a(gp);
        if ( iA != 0 ) {
            answer = answer / iA;
        }
    }
}

void
RankineMatGrad :: giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nlDamageDrivingVariable_grad, GaussPoint *gp, TimeStep *tStep)
{
    answer = nlDamageDrivingVariable_grad;
    if ( gradientDamageFormulationType == GDFT_Eikonal ) {
        double iB = this->computeEikonalInternalLength_b(gp);
        answer.times(iB);
    } else {
        answer.times(internalLength * internalLength);
    }
}

void
RankineMatGrad :: computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep)
{
    RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );
    answer =  status->giveTempCumulativePlasticStrain();
}

void
RankineMatGrad :: giveRealStressVectorGradientDamage(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep)
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

double
RankineMatGrad :: computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );
    double localCumPlastStrain = status->giveTempCumulativePlasticStrain();
    double nlCumPlastStrain = status->giveNonlocalCumulatedStrain();
    return mParam * nlCumPlastStrain + ( 1. - mParam ) * localCumPlastStrain;
}

double
RankineMatGrad :: giveNonlocalCumPlasticStrain(GaussPoint *gp)
{
    auto status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );
    return status->giveNonlocalCumulatedStrain();
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
        answer.at(1) = computeCumPlastStrain(gp, tStep);
        return 1;
    } else {
        return RankineMat :: giveIPValue(answer, gp, type, tStep);
    }
}


double
RankineMatGrad :: computeEikonalInternalLength_a(GaussPoint *gp)
{
    RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );

    double damage = status->giveTempDamage();
    return sqrt(1. - damage) * internalLength;
}

double
RankineMatGrad :: computeEikonalInternalLength_b(GaussPoint *gp)
{
    RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );

    double damage = status->giveTempDamage();
    return sqrt(1. - damage) * internalLength;
}


double
RankineMatGrad :: computeEikonalInternalLength_aPrime(GaussPoint *gp)
{
    RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );

    double damage = status->giveTempDamage();
    return -0.5 / sqrt(1. - damage) * internalLength;
}

double
RankineMatGrad :: computeEikonalInternalLength_bPrime(GaussPoint *gp)
{
    RankineMatGradStatus *status = static_cast< RankineMatGradStatus * >( this->giveStatus(gp) );

    double damage = status->giveTempDamage();
    return -0.5 / sqrt(1. - damage) * internalLength;
}

int
RankineMatGrad :: giveDimension(GaussPoint *gp)
{
    if ( gp->giveMaterialMode() == _1dMat ) {
        return 1;
    } else if ( gp->giveMaterialMode() == _PlaneStress ) {
        return 2;
    } else if ( gp->giveMaterialMode() == _PlaneStrain ) {
        return 3;
    } else if ( gp->giveMaterialMode() == _3dMat ) {
        return 3;
    } else {
        return 0;
    }
}



//=============================================================================
// GRADIENT RANKINE MATERIAL STATUS
//=============================================================================


RankineMatGradStatus :: RankineMatGradStatus(GaussPoint *g) :
    RankineMatStatus(g)
{}

void
RankineMatGradStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
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
