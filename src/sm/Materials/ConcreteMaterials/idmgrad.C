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

#include "idmgrad.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "sparsemtrx.h"
#include "Materials/isolinearelasticmaterial.h"
#include "error.h"
#include "nonlocalmaterialext.h"
#include "datastream.h"
#include "contextioerr.h"
#include "stressvector.h"
#include "strainvector.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(IsotropicGradientDamageMaterial);

IsotropicGradientDamageMaterial :: IsotropicGradientDamageMaterial(int n, Domain *d) :
    IsotropicDamageMaterial1(n, d),
    GradientDamageMaterialExtensionInterface(d)
{}


void
IsotropicGradientDamageMaterial :: initializeFrom(InputRecord &ir)
{
    IsotropicDamageMaterial1 :: initializeFrom(ir);
    GradientDamageMaterialExtensionInterface :: initializeFrom(ir);

    int formulationType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, formulationType, _IFT_IsotropicGradientDamageMaterial_formulationType);
    if ( formulationType == 0 ) {
        this->gradientDamageFormulationType = GDFT_Standard;
    } else if ( formulationType == 1 ) {
        this->gradientDamageFormulationType =   GDFT_DecreasingInteractions;
        di_rho = 0.005;
        di_eta = 5;
        IR_GIVE_OPTIONAL_FIELD(ir, di_rho, _IFT_IsotropicGradientDamageMaterial_di_rho);
        IR_GIVE_OPTIONAL_FIELD(ir, di_eta, _IFT_IsotropicGradientDamageMaterial_di_eta);
    } else if ( formulationType == 2 ) {
        this->gradientDamageFormulationType =   GDFT_Eikonal;
    } else {
        throw ValueInputException(ir, _IFT_IsotropicGradientDamageMaterial_formulationType, "Unknown gradient damage formulation");
    }


    return this->mapper.initializeFrom(ir);
}

/////////////////////////////////////////////////////////////////////////////



bool
IsotropicGradientDamageMaterial :: hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _1dMat || mode == _PlaneStress || mode == _PlaneStrain || mode == _3dMat;
}

void
IsotropicGradientDamageMaterial :: giveStiffnessMatrix(FloatMatrix &answer,
                                                       MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}


void
IsotropicGradientDamageMaterial :: giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );
    double tempDamage;
    if ( mode == ElasticStiffness ) {
        tempDamage = 0.0;
    } else {
        tempDamage = status->giveTempDamage();
        if ( tempDamage > 0.0 ) {
            tempDamage = min(tempDamage, maxOmega);
        }
    }

    this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    answer.times(1.0 - tempDamage);


    /*
     * FloatArray strainP, stressP, oldStrain, oldStress;
     * oldStress = status->giveTempStressVector();
     * oldStrain = status->giveTempStrainVector();
     * strainP = oldStrain;
     * double pert = 1.e-6 * strainP.at(1);
     * strainP.at(1) += pert;
     * double lddv;
     * double mddv = status->giveTempNonlocalDamageDrivingVariable();
     * this->giveRealStressVectorGradientDamage(stressP, lddv, gp, strainP, mddv, tStep);
     * FloatMatrix stiff;
     * stiff.resize(1,1);
     * stiff.at(1,1) = (stressP.at(1) - oldStress.at(1))/pert;
     * this->giveRealStressVectorGradientDamage(stressP, lddv, gp, oldStrain, mddv, tStep);
     */
}




void
IsotropicGradientDamageMaterial :: giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );
    FloatArray stress =  status->giveTempStressVector();
    if ( mode == TangentStiffness ) {
        double nonlocalDamageDrivingVariable =  status->giveTempNonlocalDamageDrivingVariable();
        double gPrime =  this->damageFunctionPrime(nonlocalDamageDrivingVariable, gp);
        double tempDamage = status->giveTempDamage();
        if ( tempDamage < 1. ) {
            stress.times( 1. / ( 1 - tempDamage ) );
        } else {
            stress.times(0.);
        }
        answer.initFromVector(stress, false);
        if ( tempDamage > status->giveDamage() ) {
            answer.times(-gPrime);
            /*
             * double lddv;
             * FloatArray strainP, stressP, oldStrain, oldStress;
             * oldStress = status->giveTempStressVector();
             * oldStrain = status->giveTempStrainVector();
             * double mddv = status->giveTempNonlocalDamageDrivingVariable();
             * double pert = 1.e-6 * mddv;
             * strainP = oldStrain;
             * double mddvp = mddv + pert;
             * this->giveRealStressVectorGradientDamage(stressP, lddv, gp, strainP, mddvp, tStep);
             * FloatMatrix stiff;
             * stiff.resize(1,1);
             * stiff.at(1,1) = (stressP.at(1) -oldStress.at(1))/pert;
             * this->giveRealStressVectorGradientDamage(stressP, lddv, gp, oldStrain, mddv, tStep);
             */
        } else {
            answer.times(0.);
        }
        // zero block for now
    } else {
        answer.initFromVector(stress, false);
        answer.times(0);
    }
}


void
IsotropicGradientDamageMaterial :: giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );

    FloatArray eta, reducedStrain;
    FloatArray totalStrain = status->giveTempStrainVector();
    StructuralMaterial :: giveReducedSymVectorForm( reducedStrain, totalStrain, gp->giveMaterialMode() );
    this->computeEta(eta, reducedStrain, gp, tStep);
    answer.initFromVector(eta, false);
    if ( mode == TangentStiffness ) {
        double tempDamage = status->giveTempDamage();
        if ( tempDamage > status->giveDamage() ) {
            answer.times(-1.);
            /*
             * FloatArray strainP, stressP, oldStrain, oldStress;
             * oldStress = status->giveTempStressVector();
             * oldStrain = status->giveTempStrainVector();
             * double mddv = status->giveTempNonlocalDamageDrivingVariable();
             * double lddv = status->giveTempLocalDamageDrivingVariable();
             * double lddvp;
             * strainP = oldStrain;
             * double pert = 1.e-6 * strainP.at(1);
             * strainP.at(1) += pert;
             * this->giveRealStressVectorGradientDamage(stressP, lddvp, gp, strainP, mddv, tStep);
             * FloatMatrix stiff;
             * stiff.resize(1,1);
             * stiff.at(1,1) = (lddvp - lddv) / pert;
             * this->giveRealStressVectorGradientDamage(stressP, lddv, gp, oldStrain, mddv, tStep);
             */


            if ( gradientDamageFormulationType == GDFT_Eikonal ) {
                double iA = this->computeEikonalInternalLength_a(gp);
                if ( iA != 0 ) {
                    answer.times(1. / iA);
                } else {
                    answer.times(0);
                }
            }
        }
    } else {
        answer.times(0);
    }
}



void
IsotropicGradientDamageMaterial :: giveGradientDamageStiffnessMatrix_dd_NN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    if ( gradientDamageFormulationType == GDFT_Eikonal ) {
        IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );
        double tempKappa =  status->giveTempKappa();
        double iA = this->computeEikonalInternalLength_a(gp);

        answer.resize(1, 1);
        answer.zero();

        if ( iA != 0 ) {
            answer.at(1, 1) += 1. / iA;
        }

        if ( mode == TangentStiffness ) {
            if ( tempKappa > status->giveKappa() && iA != 0 ) {
                double iAPrime = this->computeEikonalInternalLength_aPrime(gp);
                double gPrime = this->damageFunctionPrime(tempKappa, gp);
                double epsEqLocal = status->giveTempLocalDamageDrivingVariable();
                double epsEqNonLocal = status->giveTempNonlocalDamageDrivingVariable();
                answer.at(1, 1) += iAPrime / iA / iA * gPrime * ( epsEqNonLocal - epsEqLocal );
            }
        }
    } else {
        answer.clear();
    }
}

double
IsotropicGradientDamageMaterial :: computeEikonalInternalLength_a(GaussPoint *gp)
{
    IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );

    double damage = status->giveTempDamage();
    return sqrt(1. - damage) * internalLength;
}

double
IsotropicGradientDamageMaterial :: computeEikonalInternalLength_b(GaussPoint *gp)
{
    IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );

    double damage = status->giveTempDamage();
    return sqrt(1. - damage) * internalLength;
}


double
IsotropicGradientDamageMaterial :: computeEikonalInternalLength_aPrime(GaussPoint *gp)
{
    IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );

    double damage = status->giveTempDamage();
    return -0.5 / sqrt(1. - damage) * internalLength;
}

double
IsotropicGradientDamageMaterial :: computeEikonalInternalLength_bPrime(GaussPoint *gp)
{
    IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );

    double damage = status->giveTempDamage();
    return -0.5 / sqrt(1. - damage) * internalLength;
}



void
IsotropicGradientDamageMaterial :: giveGradientDamageStiffnessMatrix_dd_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    int n = this->giveDimension(gp);
    answer.resize(n, n);
    answer.beUnitMatrix();

    if ( gradientDamageFormulationType == GDFT_Standard ) {
        answer.times(internalLength * internalLength);
    } else if ( gradientDamageFormulationType == GDFT_DecreasingInteractions ) {
        double iL = this->computeInternalLength(gp);
        answer.times(iL * iL);
    } else if ( gradientDamageFormulationType == GDFT_Eikonal ) {
        double iB = this->computeEikonalInternalLength_b(gp);
        answer.times(iB);
    } else {
        OOFEM_WARNING("Unknown internalLengthDependenceType");
    }
}


void
IsotropicGradientDamageMaterial :: giveGradientDamageStiffnessMatrix_dd_BN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    if ( mode == TangentStiffness ) {
        if ( gradientDamageFormulationType == GDFT_Standard ) {
            answer.clear();
        } else if ( gradientDamageFormulationType == GDFT_DecreasingInteractions ) {
            IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );
            double damage = status->giveTempDamage();
            double tempKappa =  status->giveTempKappa();
            double iL = computeInternalLength(gp);
            double nom = -di_eta * ( 1. - di_rho ) * exp(-di_eta * damage);
            double denom = sqrt( ( 1. - exp(-di_eta) ) * ( ( 1. - di_rho ) * exp(-di_eta * damage) + di_rho - exp(-di_eta) ) );
            double gPrime = this->damageFunctionPrime(tempKappa, gp);
            double factor = iL * internalLength * nom / denom * gPrime;
            answer.initFromVector(status->giveTempNonlocalDamageDrivingVariableGrad(), false);
            if ( tempKappa > status->giveKappa() ) {
                answer.times(factor);
            } else {
                answer.times(0.);
            }
        } else if ( gradientDamageFormulationType == GDFT_Eikonal ) {
            IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );
            double tempKappa =  status->giveTempKappa();

            double iBPrime = this->computeEikonalInternalLength_bPrime(gp);
            double gPrime = this->damageFunctionPrime(tempKappa, gp);

            answer.initFromVector(status->giveTempNonlocalDamageDrivingVariableGrad(), false);

            if ( tempKappa > status->giveKappa() ) {
                answer.times(iBPrime * gPrime);
            } else {
                answer.times(0.);
            }
        } else {
            OOFEM_WARNING("Unknown internalLengthDependenceType");
        }
    } else {
        answer.clear();
    }
}


int
IsotropicGradientDamageMaterial :: giveDimension(GaussPoint *gp)
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

void
IsotropicGradientDamageMaterial :: giveNonlocalInternalForces_N_factor(double &answer, double nlDamageDrivingVariable, GaussPoint *gp, TimeStep *tStep)
{
    IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );
    answer = nlDamageDrivingVariable - status->giveTempLocalDamageDrivingVariable();

    if ( gradientDamageFormulationType == GDFT_Eikonal ) {
        double iA = this->computeEikonalInternalLength_a(gp);
        if ( iA != 0 ) {
            answer = answer / iA;
        }
    }
}

void
IsotropicGradientDamageMaterial :: giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nlDamageDrivingVariable_grad, GaussPoint *gp, TimeStep *tStep)
{
    answer = nlDamageDrivingVariable_grad;
    if ( gradientDamageFormulationType == GDFT_Eikonal ) {
        double iB = this->computeEikonalInternalLength_b(gp);
        answer.times(iB);
    } else {
        double iL = computeInternalLength(gp);
        answer.times(iL * iL);
    }
}


double
IsotropicGradientDamageMaterial :: computeInternalLength(GaussPoint *gp)
{
    if ( gradientDamageFormulationType == GDFT_Standard ) {
        return internalLength;
    } else if ( gradientDamageFormulationType == GDFT_DecreasingInteractions ) {
        IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );
        double damage = status->giveTempDamage();
        double factor = ( ( 1. - di_rho ) * exp(-di_eta * damage) + di_rho - exp(-di_eta) ) / ( 1. - exp(-di_eta) );
        return internalLength * sqrt(factor);
    } else {
        OOFEM_WARNING("Unknown internalLengthDependenceType");
        return 0.;
    }
}


void
IsotropicGradientDamageMaterial :: giveRealStressVectorGradientDamage(FloatArray &stress, double &localDamageDrivingVariable, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalDamageDrivingVariable, TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    auto status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );
    auto lmat = this->giveLinearElasticMaterial();
    FloatArray reducedTotalStrainVector, strain;

    FloatMatrix de;
    double f, damage, tempKappa = 0.0;

    this->initTempStatus(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, tStep, VM_Total);


    // compute equivalent strain
    localDamageDrivingVariable = this->computeEquivalentStrain(reducedTotalStrainVector, gp, tStep);

    if ( llcriteria == idm_strainLevelCR ) {
        // compute value of loading function if strainLevel crit apply
        f = nonlocalDamageDrivingVariable - status->giveKappa();
        if ( f <= 0.0 ) {
            // damage does not grow
            tempKappa = status->giveKappa();
            damage    = status->giveDamage();
        } else {
            // damage grow
            tempKappa = nonlocalDamageDrivingVariable;
            this->initDamaged(nonlocalDamageDrivingVariable, reducedTotalStrainVector, gp);
            // evaluate damage parameter
            damage = this->computeDamageParam(nonlocalDamageDrivingVariable, reducedTotalStrainVector, gp);
        }
    } else if ( llcriteria == idm_damageLevelCR ) {
        // evaluate damage parameter first
        tempKappa = nonlocalDamageDrivingVariable;
        this->initDamaged(nonlocalDamageDrivingVariable, strain, gp);
        damage = this->computeDamageParam(nonlocalDamageDrivingVariable, reducedTotalStrainVector, gp);
        if ( damage < status->giveDamage() ) {
            // unloading takes place
            damage = status->giveDamage();
        }
    } else {
        OOFEM_ERROR("unsupported loading/uloading criteria");
    }


    lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
    de.times(1.0 - damage);
    stress.beProductOf(de, reducedTotalStrainVector);

    // update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(stress);
    status->setTempLocalDamageDrivingVariable(localDamageDrivingVariable);
    status->setTempNonlocalDamageDrivingVariable(nonlocalDamageDrivingVariable);
    status->setTempKappa(tempKappa);
    status->setTempDamage(damage);
#ifdef keep_track_of_dissipated_energy
    status->computeWork(gp);
#endif
}




MaterialStatus *
IsotropicGradientDamageMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new IsotropicGradientDamageMaterialStatus(gp);
}


IsotropicGradientDamageMaterialStatus :: IsotropicGradientDamageMaterialStatus(GaussPoint *g) : IsotropicDamageMaterial1Status(g)
{ }


void
IsotropicGradientDamageMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    IsotropicDamageMaterial1Status :: initTempStatus();
    GradientDamageMaterialStatusExtensionInterface :: initTempStatus();
    this->tempDamage = this->damage;
}



void
IsotropicGradientDamageMaterialStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    IsotropicDamageMaterial1Status :: updateYourself(tStep);
    GradientDamageMaterialStatusExtensionInterface :: updateYourself(tStep);
}
}     // end namespace oofem
