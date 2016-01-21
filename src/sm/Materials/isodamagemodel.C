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

#include "isodamagemodel.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "dynamicinputrecord.h"

namespace oofem {
IsotropicDamageMaterial :: IsotropicDamageMaterial(int n, Domain *d) : StructuralMaterial(n, d)
    //
    // constructor
    //
{
    linearElasticMaterial = NULL;
    llcriteria = idm_strainLevelCR;
    maxOmega = 0.999999;
}


IsotropicDamageMaterial :: ~IsotropicDamageMaterial()
//
// destructor
//
{
    delete linearElasticMaterial;
}

int
IsotropicDamageMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    return mode == _3dMat || mode == _PlaneStress || mode == _PlaneStrain || mode == _1dMat;
}


void
IsotropicDamageMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                         MatResponseMode mode,
                                                         GaussPoint *gp,
                                                         TimeStep *tStep)
//
// computes full constitutive matrix for case of gp stress-strain state.
//
{
    IsotropicDamageMaterialStatus *status = static_cast< IsotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    double tempDamage;
    if ( mode == ElasticStiffness ) {
        tempDamage = 0.0;
    } else {
        tempDamage = status->giveTempDamage();
        tempDamage = min(tempDamage, maxOmega);
    }

    this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
    answer.times(1.0 - tempDamage);
    //TODO - correction for tangent mode
}


void
IsotropicDamageMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                                const FloatArray &totalStrain,
                                                TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    IsotropicDamageMaterialStatus *status = static_cast< IsotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    //StructuralCrossSection *crossSection = (StructuralCrossSection*) gp -> giveElement()->giveCrossSection();
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    FloatArray reducedTotalStrainVector;
    FloatMatrix de;
    double f, equivStrain, tempKappa = 0.0, omega = 0.0;

    this->initTempStatus(gp);
    
    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, tStep, VM_Total);

    //crossSection->giveFullCharacteristicVector(totalStrainVector, gp, reducedTotalStrainVector);

    // compute equivalent strain
    this->computeEquivalentStrain(equivStrain, reducedTotalStrainVector, gp, tStep);

    if ( llcriteria == idm_strainLevelCR ) {
        // compute value of loading function if strainLevel crit apply
        f = equivStrain - status->giveKappa();

        if ( f <= 0.0 ) {
            // damage does not grow
            tempKappa = status->giveKappa();
            omega     = status->giveDamage();
        } else {
            // damage grow
            tempKappa = equivStrain;
            this->initDamaged(tempKappa, reducedTotalStrainVector, gp);
            // evaluate damage parameter
            this->computeDamageParam(omega, tempKappa, reducedTotalStrainVector, gp);
        }
    } else if ( llcriteria == idm_damageLevelCR ) {
        // evaluate damage parameter first
        tempKappa = equivStrain;
        this->initDamaged(tempKappa, reducedTotalStrainVector, gp);
        this->computeDamageParam(omega, tempKappa, reducedTotalStrainVector, gp);
        if ( omega < status->giveDamage() ) {
            // unloading takes place
            omega = status->giveDamage();
            //printf (".");
        }
    } else {
        OOFEM_ERROR("unsupported loading/uloading criteria");
    }


    lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
    //mj
    // damage deactivation in compression for 1D model
    if ( ( reducedTotalStrainVector.giveSize() > 1 ) || ( reducedTotalStrainVector.at(1) > 0. ) ) {
        //emj
        de.times(1.0 - omega);
    }

    answer.beProductOf(de, reducedTotalStrainVector);

    // update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);
#ifdef keep_track_of_dissipated_energy
    status->computeWork(gp);
#endif
}


void IsotropicDamageMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode,
                                                         GaussPoint *gp, TimeStep *tStep)
{
    IsotropicDamageMaterialStatus *status = static_cast< IsotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    double tempDamage;
    if ( mode == ElasticStiffness ) {
        tempDamage = 0.0;
    } else {
        tempDamage = status->giveTempDamage();
        tempDamage = min(tempDamage, maxOmega);
    }

    this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    answer.times(1.0 - tempDamage);
    if ( mode == TangentStiffness ) {
        double damage = status->giveDamage();
        if ( tempDamage > damage ) {
            double tempKappa;
            FloatArray stress, strain, eta;
            FloatMatrix correctionTerm;
            stress = status->giveTempStressVector();
            strain = status->giveTempStrainVector();
            tempKappa = status->giveTempKappa();
            // effective stress
            stress.times( 1. / ( 1 - tempDamage ) );
            //Computes derivative of the equivalent strain with regards to strain
            this->computeEta(eta, strain, gp, tStep);
            //compute derivative of damage function
            double damagePrime = damageFunctionPrime(tempKappa, gp);
            // dyadic product of eff stress and eta
            correctionTerm.beDyadicProductOf(stress, eta);
            // times minus derivative of damage function
            correctionTerm.times(-damagePrime);
            // add to secant stiffness
            answer.add(correctionTerm);
        }
    }
}


void IsotropicDamageMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mode,
                                                         GaussPoint *gp, TimeStep *tStep)
{
    IsotropicDamageMaterialStatus *status = static_cast< IsotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    double tempDamage;
    if ( mode == ElasticStiffness ) {
        tempDamage = 0.0;
    } else {
        tempDamage = status->giveTempDamage();
        tempDamage = min(tempDamage, maxOmega);
    }

    this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    answer.times(1.0 - tempDamage);
    //TODO - correction for tangent mode
}


void IsotropicDamageMaterial :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode,
                                                      GaussPoint *gp, TimeStep *tStep)
{
    IsotropicDamageMaterialStatus *status = static_cast< IsotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    double tempDamage;
    if ( mode == ElasticStiffness ) {
        tempDamage = 0.0;
    } else {
        tempDamage = status->giveTempDamage();
        tempDamage = min(tempDamage, maxOmega);
    }

    this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    answer.times(1.0 - tempDamage);
    //TODO - correction for tangent mode
}

#ifdef __OOFEG
#endif


int
IsotropicDamageMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    IsotropicDamageMaterialStatus *status = static_cast< IsotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveDamage();
        return 1;
    } else if ( type == IST_DamageTensor ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveDamage();
        return 1;
    } else if ( type == IST_PrincipalDamageTensor ) {
        answer.resize(3);
        answer.zero();
        answer.at(1) = status->giveDamage();
        return 1;
    } else if ( type == IST_DamageTensorTemp ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveTempDamage();
        return 1;
    } else if ( type == IST_PrincipalDamageTempTensor ) {
        answer.resize(3);
        answer.zero();
        answer.at(1) = status->giveTempDamage();
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveKappa();
        return 1;
    } else if ( type == IST_CharacteristicLength ) {
        answer.resize(1);
        answer.at(1) = status->giveLe();
        return 1;
    } else if ( type == IST_CrackDirs ) {
        answer.resize(1);
        answer.at(1) = status->giveCrackAngle();
        return 1;
    } else if ( type == IST_CrackVector ) {
        status->giveCrackVector(answer);
        return 1;

#ifdef keep_track_of_dissipated_energy
    } else if ( type == IST_StressWorkDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveStressWork();
        return 1;
    } else if ( type == IST_DissWorkDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveDissWork();
    } else if ( type == IST_FreeEnergyDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveStressWork() - status->giveDissWork();
        return 1;

#endif
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }

    return 1; // to make the compiler happy
}


void
IsotropicDamageMaterial :: giveThermalDilatationVector(FloatArray &answer,
                                                       GaussPoint *gp,  TimeStep *tStep)
//
// returns a FloatArray(6) of initial strain vector
// eps_0 = {exx_0, eyy_0, ezz_0, gyz_0, gxz_0, gxy_0}^T
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    answer.resize(6);
    answer.zero();
    answer.at(1) = this->tempDillatCoeff;
    answer.at(2) = this->tempDillatCoeff;
    answer.at(3) = this->tempDillatCoeff;
}

double IsotropicDamageMaterial :: give(int aProperty, GaussPoint *gp)
{
    return linearElasticMaterial->give(aProperty, gp);
}

IRResultType
IsotropicDamageMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    //Set limit on the maximum isotropic damage parameter if needed
    IR_GIVE_OPTIONAL_FIELD(ir, maxOmega, _IFT_IsotropicDamageMaterial_maxOmega);
    maxOmega = min(maxOmega, 0.999999);
    maxOmega = max(maxOmega, 0.0);

    IR_GIVE_FIELD(ir, tempDillatCoeff, _IFT_IsotropicDamageMaterial_talpha);
    return StructuralMaterial :: initializeFrom(ir);
}


void
IsotropicDamageMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);
    input.setField(this->maxOmega, _IFT_IsotropicDamageMaterial_maxOmega);
    input.setField(this->tempDillatCoeff, _IFT_IsotropicDamageMaterial_talpha);
}



IsotropicDamageMaterialStatus :: IsotropicDamageMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    kappa = tempKappa = 0.0;
    damage = tempDamage = 0.0;
    le = 0.0;
    crack_angle = -1000.0;
    crackVector.resize(3);
    crackVector.zero();
#ifdef keep_track_of_dissipated_energy
    stressWork = tempStressWork = 0.0;
    dissWork = tempDissWork = 0.0;
#endif
}


IsotropicDamageMaterialStatus :: ~IsotropicDamageMaterialStatus()
{ }


void
IsotropicDamageMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( this->kappa > 0 && this->damage <= 0 ) {
        fprintf(file, "kappa %f", this->kappa);
    } else if ( this->damage > 0.0 ) {
        fprintf( file, "kappa %f, damage %f crackVector %f %f %f", this->kappa, this->damage, this->crackVector.at(1), this->crackVector.at(2), this->crackVector.at(3) );

#ifdef keep_track_of_dissipated_energy
        fprintf(file, ", dissW %f, freeE %f, stressW %f ", this->dissWork, ( this->stressWork ) - ( this->dissWork ), this->stressWork);
    } else {
        fprintf(file, "stressW %f ", this->stressWork);
#endif
    }

    fprintf(file, "}\n");
}


void
IsotropicDamageMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
    //mj 14 July 2010 - should be discussed with Borek !!!
    //this->tempDamage = this->damage;
#ifdef keep_track_of_dissipated_energy
    this->tempStressWork = this->stressWork;
    this->tempDissWork = this->dissWork;
#endif
}


void
IsotropicDamageMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    this->kappa = this->tempKappa;
    this->damage = this->tempDamage;
#ifdef keep_track_of_dissipated_energy
    this->stressWork = this->tempStressWork;
    this->dissWork = this->tempDissWork;
#endif
}

void
IsotropicDamageMaterialStatus :: giveCrackVector(FloatArray &answer)
{
    answer = crackVector;
    answer.times(damage);
}


contextIOResultType
IsotropicDamageMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data
    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream.write(stressWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(dissWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#endif

    return CIO_OK;
}

contextIOResultType
IsotropicDamageMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream.read(stressWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(dissWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#endif

    return CIO_OK;
}

#ifdef keep_track_of_dissipated_energy
void
IsotropicDamageMaterialStatus :: computeWork(GaussPoint *gp)
{
    // strain increment
    FloatArray deps;
    deps.beDifferenceOf(tempStrainVector, strainVector);

    // increment of stress work density
    double dSW = ( tempStressVector.dotProduct(deps) + stressVector.dotProduct(deps) ) / 2.;
    tempStressWork = stressWork + dSW;

    // elastically stored energy density
    double We = tempStressVector.dotProduct(tempStrainVector) / 2.;

    // dissipative work density
    tempDissWork = tempStressWork - We;
}
#endif
} // end namespace oofem
