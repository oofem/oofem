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

#include "isodamagemodel.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"

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
    if ( ( mode == _3dMat ) || ( mode == _PlaneStress ) ||
        ( mode == _PlaneStrain ) || ( mode == _1dMat ) ) {
        return 1;
    }

    return 0;
}


void
IsotropicDamageMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                         MatResponseForm form,
                                                         MatResponseMode mode,
                                                         GaussPoint *gp,
                                                         TimeStep *atTime)
//
// computes full constitutive matrix for case of gp stress-strain state.
//
{
    IsotropicDamageMaterialStatus *status = ( IsotropicDamageMaterialStatus * ) this->giveStatus(gp);
    double om;
    if ( mode == ElasticStiffness ) {
        om = 0.0;
    } else {
        om = status->giveTempDamage();
        om = min(om, maxOmega);
    }

    this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
    answer.times(1.0 - om);
}


void
IsotropicDamageMaterial :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                                const FloatArray &totalStrain,
                                                TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    IsotropicDamageMaterialStatus *status = ( IsotropicDamageMaterialStatus * ) this->giveStatus(gp);
    //StructuralCrossSection *crossSection = (StructuralCrossSection*) gp -> giveElement()->giveCrossSection();
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    FloatArray strainVector, reducedTotalStrainVector;
    FloatMatrix de;
    double f, equivStrain, tempKappa = 0.0, omega = 0.0;

    this->initGpForNewStep(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, atTime, VM_Total);

    //crossSection->giveFullCharacteristicVector(totalStrainVector, gp, reducedTotalStrainVector);

    // compute equivalent strain
    this->computeEquivalentStrain(equivStrain, reducedTotalStrainVector, gp, atTime);

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
        _error("giveRealStressVector: unsupported loading/uloading criteria");
    }


    lmat->giveCharacteristicMatrix(de, ReducedForm, SecantStiffness, gp, atTime);
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


void IsotropicDamageMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode,
                                                         GaussPoint *gp, TimeStep *atTime)
{
    IsotropicDamageMaterialStatus *status = ( IsotropicDamageMaterialStatus * ) this->giveStatus(gp);
    double om;
    if ( mode == ElasticStiffness ) {
        om = 0.0;
    } else {
        om = status->giveTempDamage();
        om = min(om, maxOmega);
    }

    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    answer.times(1.0 - om);
}


void IsotropicDamageMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode,
                                                         GaussPoint *gp, TimeStep *atTime)
{
    IsotropicDamageMaterialStatus *status = ( IsotropicDamageMaterialStatus * ) this->giveStatus(gp);
    double om;
    if ( mode == ElasticStiffness ) {
        om = 0.0;
    } else {
        om = status->giveTempDamage();
        om = min(om, maxOmega);
    }

    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    answer.times(1.0 - om);
}


void IsotropicDamageMaterial :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode,
                                                      GaussPoint *gp, TimeStep *atTime)
{
    IsotropicDamageMaterialStatus *status = ( IsotropicDamageMaterialStatus * ) this->giveStatus(gp);
    double om;
    if ( mode == ElasticStiffness ) {
        om = 0.0;
    } else {
        om = status->giveTempDamage();
        om = min(om, maxOmega);
    }

    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    answer.times(1.0 - om);
}

#ifdef __OOFEG
#endif


int
IsotropicDamageMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    IsotropicDamageMaterialStatus *status = ( IsotropicDamageMaterialStatus * ) this->giveStatus(aGaussPoint);
    if ( ( type == IST_DamageTensor ) || ( type == IST_PrincipalDamageTensor ) ) {
        answer.resize(1);
        answer.at(1) = status->giveDamage();
        return 1;
    } else if ( ( type == IST_DamageTensorTemp ) || ( type == IST_PrincipalDamageTempTensor ) ) {
        answer.resize(1);
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
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }

    return 1; // to make the compiler happy
}


InternalStateValueType
IsotropicDamageMaterial :: giveIPValueType(InternalStateType type)
{
    if ( ( type == IST_DamageTensor ) || ( type == IST_DamageTensorTemp ) ||
        ( type == IST_PrincipalDamageTensor ) || ( type == IST_PrincipalDamageTempTensor ) ) {
        return ISVT_TENSOR_S3;
    } else if ( type == IST_MaxEquivalentStrainLevel || type == IST_CharacteristicLength || type == IST_CrackDirs ) {
        return ISVT_SCALAR;

#ifdef keep_track_of_dissipated_energy
    } else if ( type == IST_DissWorkDensity ) {
        return ISVT_SCALAR;
    } else if ( type == IST_StressWorkDensity ) {
        return ISVT_SCALAR;
    } else if ( type == IST_FreeEnergyDensity ) {
        return ISVT_SCALAR;

#endif
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }
}


int
IsotropicDamageMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( ( type == IST_DamageTensor ) || ( type == IST_DamageTensorTemp ) ||
        ( type == IST_PrincipalDamageTensor ) || ( type == IST_PrincipalDamageTempTensor ) ) {
        answer.resize(9);
        answer.at(1) = 1;
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel || type == IST_CharacteristicLength || type == IST_CrackDirs ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;

#ifdef keep_track_of_dissipated_energy
    } else if ( ( type == IST_DissWorkDensity ) || ( type == IST_StressWorkDensity ) || ( type == IST_FreeEnergyDensity ) ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;

#endif
    } else {
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}


int
IsotropicDamageMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( ( type == IST_DamageTensor ) || ( type == IST_DamageTensorTemp ) ||
        ( type == IST_PrincipalDamageTensor ) || ( type == IST_PrincipalDamageTempTensor ) ||
        ( type == IST_MaxEquivalentStrainLevel || type == IST_CharacteristicLength || type == IST_CrackDirs ) ) {
        return 1;

#ifdef keep_track_of_dissipated_energy
    } else if ( ( type == IST_StressWorkDensity ) ||
               ( type == IST_DissWorkDensity ) || ( type == IST_FreeEnergyDensity ) ) {
        return 1;

#endif
    } else {
        return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
    }
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


IRResultType
IsotropicDamageMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    //Set limit on the maximum isotropic damage parameter if needed
    IR_GIVE_OPTIONAL_FIELD(ir, maxOmega, IFT_IsotropicDamageMaterial_maxOmega, "maxomega"); // Macro
    maxOmega = min(maxOmega, 0.999999);
    maxOmega = max(maxOmega, 0.0);

    IR_GIVE_FIELD(ir, tempDillatCoeff, IFT_IsotropicDamageMaterial_talpha, "talpha"); // Macro
    return StructuralMaterial :: initializeFrom(ir);
}


int
IsotropicDamageMaterial :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    StructuralMaterial :: giveInputRecordString(str, keyword);
    sprintf(buff, " talpha %e", this->tempDillatCoeff);
    str += buff;

    return 1;
}



IsotropicDamageMaterialStatus :: IsotropicDamageMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    kappa = tempKappa = 0.0;
    damage = tempDamage = 0.0;
    le = 0.0;
    crack_angle = -1000.0;
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
    if ( this->damage > 0.0 ) {
        fprintf(file, "kappa %f, damage %f ", this->kappa, this->damage);

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
IsotropicDamageMaterialStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
    this->kappa = this->tempKappa;
    this->damage = this->tempDamage;
#ifdef keep_track_of_dissipated_energy
    this->stressWork = this->tempStressWork;
    this->dissWork = this->tempDissWork;
#endif
}


contextIOResultType
IsotropicDamageMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data
    if ( !stream->write(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream->write(& stressWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& dissWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#endif

    return CIO_OK;
}

contextIOResultType
IsotropicDamageMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream->read(& stressWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& dissWork, 1) ) {
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
