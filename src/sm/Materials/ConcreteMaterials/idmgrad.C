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

IsotropicGradientDamageMaterial :: IsotropicGradientDamageMaterial(int n, Domain *d) : IsotropicDamageMaterial1(n, d), GradientDamageMaterialExtensionInterface(d)
    //
    // constructor
    //
{}


IsotropicGradientDamageMaterial :: ~IsotropicGradientDamageMaterial()
//
// destructor
//
{ }

IRResultType
IsotropicGradientDamageMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    result = IsotropicDamageMaterial1 :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    result = GradientDamageMaterialExtensionInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    int internalLengthDependence = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, internalLengthDependence, _IFT_IsotropicGradientDamageMaterial_internalLengthDependence);
    if ( internalLengthDependence == 0 ) {
        this->internalLengthDependenceType = ILD_None;
    } else if ( internalLengthDependence == 1 ) {
        this->internalLengthDependenceType =  ILD_Damage_DecreasingInteractions;
        di_rho = 0.005;
        di_eta = 5;
        IR_GIVE_OPTIONAL_FIELD(ir, di_rho, _IFT_IsotropicGradientDamageMaterial_di_rho);
        IR_GIVE_OPTIONAL_FIELD(ir, di_eta, _IFT_IsotropicGradientDamageMaterial_di_eta);
    } else {
        OOFEM_WARNING("Unknown internalLenghtDependence %d", internalLengthDependence);
        return IRRT_BAD_FORMAT;
    }


    return this->mapper.initializeFrom(ir);
}

/////////////////////////////////////////////////////////////////////////////



int
IsotropicGradientDamageMaterial :: hasMaterialModeCapability(MaterialMode mode)
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
}




void
IsotropicGradientDamageMaterial :: giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );
    FloatArray stress =  status->giveTempStressVector();
    if ( mode == TangentStiffness ) {
        double tempKappa =  status->giveTempKappa();
        double gPrime =  this->damageFunctionPrime(tempKappa, gp);
        double tempDamage = status->giveTempDamage();
        if ( tempDamage < 1. ) {
            stress.times( 1. / ( 1 - tempDamage ) );
        } else {
            stress.times(0.);
        }
        answer.initFromVector(stress, false);
        if ( tempDamage > status->giveDamage() ) {
            answer.times(-gPrime);
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

    FloatArray eta;
    FloatArray totalStrain = status->giveTempStrainVector();
    this->computeEta(eta, totalStrain, gp, tStep);
    answer.initFromVector(eta, false);
    answer.times(-1.);
    if ( mode != TangentStiffness ) {
        // zero block for now
        answer.times(0);
    }
}



void
IsotropicGradientDamageMaterial :: giveGradientDamageStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(1, 1);
    answer.at(1, 1) = 1.;
}


void
IsotropicGradientDamageMaterial :: giveGradientDamageStiffnessMatrix_dd_l(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    int n = this->giveDimension(gp);
    answer.resize(n, n);
    answer.beUnitMatrix();

    if ( internalLengthDependenceType == ILD_None ) {
        answer.times(internalLength * internalLength);
    } else if ( internalLengthDependenceType == ILD_Damage_DecreasingInteractions ) {
        double iL = this->computeInternalLength(gp);
        answer.times(iL * iL);
    } else {
        OOFEM_WARNING("Unknown internalLengthDependenceType");
    }
}


void
IsotropicGradientDamageMaterial :: giveGradientDamageStiffnessMatrix_dd_dl(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    if ( mode == TangentStiffness ) {
        if ( internalLengthDependenceType == ILD_None ) {
            answer.clear();
        } else if ( internalLengthDependenceType == ILD_Damage_DecreasingInteractions ) {
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
}

void
IsotropicGradientDamageMaterial :: giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nlDamageDrivingVariable_grad, GaussPoint *gp, TimeStep *tStep)
{
    answer = nlDamageDrivingVariable_grad;
    double iL = computeInternalLength(gp);
    answer.times(iL * iL);
}


double
IsotropicGradientDamageMaterial :: computeInternalLength(GaussPoint *gp)
{
    if ( internalLengthDependenceType == ILD_None ) {
        return internalLength;
    } else if ( internalLengthDependenceType == ILD_Damage_DecreasingInteractions ) {
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
    IsotropicGradientDamageMaterialStatus *status = static_cast< IsotropicGradientDamageMaterialStatus * >( this->giveStatus(gp) );
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    FloatArray reducedTotalStrainVector, strain;

    FloatMatrix de;
    double f, damage, tempKappa = 0.0;

    this->initTempStatus(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, tStep, VM_Total);


    // compute equivalent strain
    this->computeEquivalentStrain(localDamageDrivingVariable, reducedTotalStrainVector, gp, tStep);

    if ( llcriteria == idm_strainLevelCR ) {
        // compute value of loading function if strainLevel crit apply
        f = nonlocalDamageDrivingVariable - status->giveKappa();
        if ( f <= 0.0 ) {
            // damage does not grow
            tempKappa = status->giveKappa();
            damage     = status->giveDamage();
        } else {
            // damage grow
            tempKappa = nonlocalDamageDrivingVariable;
            this->initDamaged(nonlocalDamageDrivingVariable, reducedTotalStrainVector, gp);
            // evaluate damage parameter
            this->computeDamageParam(damage, nonlocalDamageDrivingVariable, reducedTotalStrainVector, gp);
        }
    } else if ( llcriteria == idm_damageLevelCR ) {
        // evaluate damage parameter first
        tempKappa = nonlocalDamageDrivingVariable;
        this->initDamaged(nonlocalDamageDrivingVariable, strain, gp);
        this->computeDamageParam(damage, nonlocalDamageDrivingVariable, reducedTotalStrainVector, gp);
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


IsotropicGradientDamageMaterialStatus :: ~IsotropicGradientDamageMaterialStatus()
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
