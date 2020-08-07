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

#include "concretedpm2.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "intarray.h"
#include "datastream.h"
#include "contextioerr.h"
#include "timestep.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "mathfem.h"
#include "classfactory.h"
#include <limits>

namespace oofem {
REGISTER_Material(ConcreteDPM2);

ConcreteDPM2Status::ConcreteDPM2Status(GaussPoint *gp) :
    StructuralMaterialStatus(gp)
{
    strainVector.resize(6);
    tempStrainVector.resize(6);
    stressVector.resize(6);
    tempStressVector.resize(6);
}

void
ConcreteDPM2Status::initTempStatus()
{
    // Call the function of the parent class to initialize the variables defined there.
    StructuralMaterialStatus::initTempStatus();

    tempReducedStrain = reducedStrain;
    tempPlasticStrain = plasticStrain;

    tempKappaP = kappaP;
    tempAlpha = alpha;

    //Damage part
    tempEquivStrain = equivStrain;

    tempEquivStrainTension = equivStrainTension;
    tempEquivStrainCompression = equivStrainCompression;

    tempKappaDTension = kappaDTension;
    tempKappaDCompression = kappaDCompression;

    tempKappaDTensionOne = kappaDTensionOne;
    tempKappaDCompressionOne = kappaDCompressionOne;
    tempKappaDTensionTwo = kappaDTensionTwo;
    tempKappaDCompressionTwo = kappaDCompressionTwo;

    tempDamageTension = damageTension;
    tempDamageCompression = damageCompression;

    temp_state_flag = state_flag;
    tempRateFactor = rateFactor;
    tempRateStrain = rateStrain;

    tempReturnType = returnType;
    tempReturnResult = returnResult;

#ifdef keep_track_of_dissipated_energy
    tempStressWork = stressWork;
    tempDissWork = dissWork;
#endif
}

void
ConcreteDPM2Status::updateYourself(TimeStep *tStep)
{
    // Call corresponding function of the parent class to update
    // variables defined there.
    StructuralMaterialStatus::updateYourself(tStep);

    // update variables defined in ConcreteDPM2Status

    alpha = tempAlpha;

    reducedStrain = tempReducedStrain;

    //Plasticity part
    plasticStrain = tempPlasticStrain;
    kappaP = tempKappaP;


    //Damage part
    equivStrain = tempEquivStrain;
    equivStrainTension = tempEquivStrainTension;
    equivStrainCompression = tempEquivStrainCompression;

    kappaDTension = tempKappaDTension;
    kappaDCompression = tempKappaDCompression;

    kappaDTensionOne = tempKappaDTensionOne;
    kappaDCompressionOne = tempKappaDCompressionOne;

    kappaDTensionTwo = tempKappaDTensionTwo;
    kappaDCompressionTwo = tempKappaDCompressionTwo;

    damageTension = tempDamageTension;
    damageCompression = tempDamageCompression;

    state_flag = temp_state_flag;

    rateFactor = tempRateFactor;

    rateStrain = tempRateStrain;

    returnType = tempReturnType;

    returnResult = tempReturnResult;

#ifdef keep_track_of_dissipated_energy
    stressWork = tempStressWork;
    dissWork = tempDissWork;
#endif
}

void
ConcreteDPM2Status::printOutputAt(FILE *file, TimeStep *tStep) const
{
    // Call corresponding function of the parent class to print
    StructuralMaterialStatus::printOutputAt(file, tStep);

    fprintf(file, "\tstatus { ");

    // print status flag
    switch ( state_flag ) {
    case ConcreteDPM2Status::ConcreteDPM2_Elastic:
        fprintf(file, "Elastic, ");
        break;
    case ConcreteDPM2Status::ConcreteDPM2_Unloading:
        fprintf(file, "Unloading, ");
        break;
    case ConcreteDPM2Status::ConcreteDPM2_Plastic:
        fprintf(file, "Plastic, ");
        break;
    case ConcreteDPM2Status::ConcreteDPM2_Damage:
        fprintf(file, "Damage, ");
        break;
    case ConcreteDPM2Status::ConcreteDPM2_PlasticDamage:
        fprintf(file, "PlasticDamage, ");
        break;
    }

    // print plastic strain vector and inelastic strain vector
    auto inelasticStrainVector = ( ( FloatArrayF< 6 >(strainVector) - plasticStrain ) * damageTension + plasticStrain ) * le;

    fprintf(file, " reduced ");
    for ( auto &val : reducedStrain ) {
        fprintf(file, " %.10e", val);
    }

    fprintf(file, " plastic ");
    for ( auto &val : plasticStrain ) {
        fprintf(file, " %.10e", val);
    }

    fprintf(file, " inelastic ");
    for ( auto &val : inelasticStrainVector ) {
        fprintf(file, " %.10e", val);
    }

    fprintf(file, " equivStrain %.10e,", equivStrain);

    fprintf(file, " kappaDTension %.10e,", kappaDTension);

    fprintf(file, " kappaDCompression %.10e,", kappaDCompression);

    fprintf(file, " kappaP %.10e,", kappaP);

    fprintf(file, " kappaDTensionOne %.10e,", kappaDTensionOne);

    fprintf(file, " kappaDCompressionOne %.10e,", kappaDCompressionOne);

    fprintf(file, " kappaDTensionTwo %.10e,", kappaDTensionTwo);

    fprintf(file, " kappaDCompressionTwo %.10e,", kappaDCompressionTwo);

    fprintf(file, " damageTension %.10e,", damageTension);

    fprintf(file, " damageCompression %.10e,", damageCompression);

    fprintf(file, " alpha %.10e,", this->alpha);

#ifdef keep_track_of_dissipated_energy
    fprintf(file, " dissW %g, freeE %g, stressW %g ", this->dissWork, ( this->stressWork ) - ( this->dissWork ), this->stressWork);
#endif
    fprintf(file, "}\n");
}

void
ConcreteDPM2Status::saveContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus::saveContext(stream, mode);

    contextIOResultType iores;


    if ( ( iores = plasticStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    if ( ( iores = reducedStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( !stream.write(kappaP) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    if ( !stream.write(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(alpha) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(equivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(equivStrainTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(equivStrainCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDTensionOne) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDCompressionOne) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDTensionTwo) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDCompressionTwo) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    if ( !stream.write(damageTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(damageCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(deltaEquivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(rateFactor) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(rateStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(state_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(returnType) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(returnResult) ) {
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
}


void
ConcreteDPM2Status::restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus::restoreContext(stream, mode);

    contextIOResultType iores;

    if ( ( iores = plasticStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    if ( ( iores = reducedStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(kappaP) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(alpha) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(equivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(equivStrainTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(equivStrainCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappaDTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappaDCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappaDTensionOne) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappaDCompressionOne) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappaDTensionTwo) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappaDCompressionTwo) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    if ( !stream.read(damageTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(damageCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(deltaEquivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(rateFactor) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(rateStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(state_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(returnType) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(returnResult) ) {
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
}

#ifdef keep_track_of_dissipated_energy
void
ConcreteDPM2Status::computeWork(GaussPoint *gp, double gf)
{
    auto tempTotalstrain = tempReducedStrain;
    auto totalstrain = reducedStrain;

    //Calculate increase or decrease of total strain tensor during iteration/step
    auto deltaTotalStrain = tempTotalstrain - totalstrain;

    //Ask for stress tensor at step
    auto stress = tempStressVector;

    //Calculate increase/decrease in total work
    double dSW = ( dot(tempStressVector, deltaTotalStrain) + dot(stressVector, deltaTotalStrain) ) / 2.;

    double tempStressWork = this->giveTempStressWork() + dSW;

    //Calculate temporary elastic strain
    auto tempElasticStrain = tempTotalstrain - tempPlasticStrain;

    //Calculate elastically stored energy density
    double We = dot(tempStressVector, tempElasticStrain) / 2.;

    // dissipative work density
    tempDissWork = tempStressWork - We;

    // to avoid extremely small negative dissipation due to round-off error
    // (note: gf is the dissipation density at complete failure, per unit volume)
    //Rough estimation of gf which is ok for this purpose

    if ( fabs(tempDissWork) < 1.e-12 * gf ) {
        tempDissWork = 0.;
    }

    this->setTempStressWork(tempStressWork);
    this->setTempDissWork(tempDissWork);
}
#endif

//   ********************************
//   *** CLASS CONCRETE DAMAGE PLASTICITY MODEL 2 ***
//   ********************************

#define IDM_ITERATION_LIMIT 1.e-8

ConcreteDPM2::ConcreteDPM2(int n, Domain *d) :
    StructuralMaterial(n, d),
    linearElasticMaterial(n, d)
{}


bool
ConcreteDPM2::hasMaterialModeCapability(MaterialMode mode) const
//
// returns whether receiver supports given mode
//
{
    return mode == _3dMat;
}


void
ConcreteDPM2::initializeFrom(InputRecord &ir)
{
    // call the corresponding service for the linear elastic material
    StructuralMaterial::initializeFrom(ir);
    linearElasticMaterial.initializeFrom(ir);

    //damage flag
    this->damageFlag = 1; //Default value using damage in tension and compression according to IJSS CDPM2 paper.
    IR_GIVE_OPTIONAL_FIELD(ir, this->damageFlag, _IFT_ConcreteDPM2_damflag);

    // elastic parameters
    IR_GIVE_FIELD(ir, this->eM, _IFT_IsotropicLinearElasticMaterial_e)
    IR_GIVE_FIELD(ir, this->nu, _IFT_IsotropicLinearElasticMaterial_n);
    propertyDictionary.add('E', this->eM);
    propertyDictionary.add('n', this->nu);

    this->gM = this->eM / ( 2. * ( 1. + this->nu ) );
    this->kM = this->eM / ( 3. * ( 1. - 2. * this->nu ) );

    IR_GIVE_FIELD(ir, this->fc, _IFT_ConcreteDPM2_fc);
    IR_GIVE_FIELD(ir, this->ft, _IFT_ConcreteDPM2_ft);

    this->e0 = this->ft / this->eM;

    // default parameters
    this->ecc = 0.525;
    IR_GIVE_OPTIONAL_FIELD(ir, this->ecc, _IFT_ConcreteDPM2_ecc);
    this->yieldHardInitial = 0.3;
    IR_GIVE_OPTIONAL_FIELD(ir, this->yieldHardInitial, _IFT_ConcreteDPM2_kinit);

    //Inclination at transition point
    this->yieldHardPrimePeak = 0.5;
    IR_GIVE_OPTIONAL_FIELD(ir, this->yieldHardPrimePeak, _IFT_ConcreteDPM2_hp);

    if ( this->yieldHardPrimePeak < 0 ) {
        this->yieldHardPrimePeak = 0.;
        OOFEM_WARNING("kPrimePeak cannot be less than zero");
    } else if ( this->yieldHardPrimePeak > ( 1. - this->yieldHardInitial ) ) {
        this->yieldHardPrimePeak = 1. - this->yieldHardInitial;
        OOFEM_WARNING("kPrimePeak cannot be greater than 1.-kinit");
    }

    this->AHard = 8.e-2;
    IR_GIVE_OPTIONAL_FIELD(ir, this->AHard, _IFT_ConcreteDPM2_ahard);
    this->BHard = 3.e-3;
    IR_GIVE_OPTIONAL_FIELD(ir, this->BHard, _IFT_ConcreteDPM2_bhard);
    this->CHard = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->CHard, _IFT_ConcreteDPM2_chard);
    this->DHard = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->DHard, _IFT_ConcreteDPM2_dhard);
    this->dilationConst = 0.85;
    IR_GIVE_OPTIONAL_FIELD(ir, this->dilationConst, _IFT_ConcreteDPM2_dilation);

    this->softeningType = 1; //0-Linear softening; 1-Bilinear softening; 2-Exponential
    IR_GIVE_OPTIONAL_FIELD(ir, this->softeningType, _IFT_ConcreteDPM2_softeningType);

    if ( softeningType > 2 ) {
        throw ValueInputException(ir, _IFT_ConcreteDPM2_softeningType, "softening type not implemented");
    }

    IR_GIVE_FIELD(ir, this->wf, _IFT_ConcreteDPM2_wf);

    if ( this->softeningType == 1 ) {
        this->ftOne = 0.3 * this->ft;
        IR_GIVE_OPTIONAL_FIELD(ir, this->ftOne, _IFT_ConcreteDPM2_ftOne);
        this->wfOne = 0.15 * this->wf;
        IR_GIVE_OPTIONAL_FIELD(ir, this->wfOne, _IFT_ConcreteDPM2_wfOne);
    }

    this->efCompression = 100.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->efCompression, _IFT_ConcreteDPM2_efc);

    this->ASoft = 15;
    IR_GIVE_OPTIONAL_FIELD(ir, this->ASoft, _IFT_ConcreteDPM2_asoft);

    this->helem = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->helem, _IFT_ConcreteDPM2_helem);


    //Compute m
    m = 3. * ( pow(this->fc, 2.) - pow(this->ft, 2.) ) / ( this->fc * this->ft ) * this->ecc / ( this->ecc + 1. );

    //Compute default value of dilationConst
    this->yieldTol = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->yieldTol, _IFT_ConcreteDPM2_yieldtol);

    this->yieldTolDamage = this->yieldTol * 10.;

    this->newtonIter = 100;
    IR_GIVE_OPTIONAL_FIELD(ir, this->newtonIter, _IFT_ConcreteDPM2_newtoniter);

    this->strengthRateType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->strengthRateType, _IFT_ConcreteDPM2_strengthratetype);
    if ( this->strengthRateType < 0 || this->strengthRateType > 2 ) {
        OOFEM_ERROR("strength rate type not implemented. Must be 0, 1 or 2\n");
    }

    this->deltaTime = -1.;
    this->energyRateType = 0;
    if ( this->strengthRateType > 0 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, this->deltaTime, _IFT_ConcreteDPM2_deltatime);
        IR_GIVE_OPTIONAL_FIELD(ir, this->energyRateType, _IFT_ConcreteDPM2_energyratetype);
    }
}


void
ConcreteDPM2::saveContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterial::saveContext(stream, mode);
    if ( ( mode & CM_Definition ) ) {
        if ( !stream.write(fc) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(ft) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(ecc) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(damageFlag) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(e0) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(AHard) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(BHard) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(CHard) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(DHard) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(hardeningModulus) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(ASoft) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(yieldHardPrimePeak) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(yieldHardInitial) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(dilationConst) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(m) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(mQ) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(helem) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(eM) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(gM) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(kM) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(nu) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(efCompression) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(wf) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(wfOne) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(ftOne) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(yieldTol) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(yieldTolDamage) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(newtonIter) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(softeningType) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(deltaTime) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(strengthRateType) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(energyRateType) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }


    linearElasticMaterial.saveContext(stream, mode);
}

void
ConcreteDPM2::restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterial::restoreContext(stream, mode);
    if ( ( mode & CM_Definition ) ) {
        if ( !stream.read(fc) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(ft) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(ecc) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(damageFlag) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(e0) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(AHard) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(BHard) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(CHard) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(DHard) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(hardeningModulus) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(ASoft) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(yieldHardPrimePeak) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(yieldHardInitial) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(dilationConst) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(m) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(mQ) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(helem) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(eM) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(gM) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(kM) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(nu) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(efCompression) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(wf) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(wfOne) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(ftOne) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(yieldTol) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(yieldTolDamage) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(newtonIter) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(softeningType) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(deltaTime) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(strengthRateType) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(energyRateType) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
    linearElasticMaterial.restoreContext(stream, mode);
}




FloatArrayF< 6 >
ConcreteDPM2::giveRealStressVector_3d(const FloatArrayF< 6 > &fullStrainVector, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    // Initialize temp variables for this gauss point
    status->initTempStatus();

    //    Remove thermal/shrinkage strains
    auto thermalStrain = this->computeStressIndependentStrainVector_3d(gp, tStep, VM_Total);

    auto strainVector = fullStrainVector - thermalStrain;

    status->letTempReducedStrainBe(strainVector);

    //Calculate time increment
    double dt = deltaTime;
    if ( dt == -1 ) {
        if ( tStep->giveTimeIncrement() == 0 ) { //Problem with the first step. For some reason the time increment is zero
            dt = 1.;
        } else {
            dt = tStep->giveTimeIncrement();
        }
    }

    auto D = this->linearElasticMaterial.give3dMaterialStiffnessMatrix(ElasticStiffness, gp, tStep);
    // perform plasticity return
    auto effectiveStress = performPlasticityReturn(gp, D, strainVector);

    FloatArrayF< 6 >effectiveStressTension;
    FloatArrayF< 6 >effectiveStressCompression;
    FloatArrayF< 6 >stress;

    double alpha = 0.;

    if ( this->damageFlag != 0 ) {//Apply damage
        alpha  = computeAlpha(effectiveStressTension, effectiveStressCompression, effectiveStress);
        auto damages = computeDamage(strainVector, D, dt, gp, tStep, alpha, effectiveStress);

        if ( this->damageFlag == 1 ) { //Default as described in IJSS CDPM2 article
            stress = effectiveStressTension * ( 1. - damages.at(1) ) + effectiveStressCompression * ( 1. - damages.at(2) );
        } else if ( this->damageFlag == 2 ) {  //Simplified version without split of stress but two damage variables
            stress = effectiveStress * ( 1. - ( 1. - alpha ) * damages.at(1) ) * ( 1. - alpha * damages.at(2) );
        } else if ( this->damageFlag == 3 ) { //Consider only tensile damage. Reduction to a fully isotropic model. Similar to CDPM article.
            stress = effectiveStress * ( 1. - damages.at(1) );
        } else {
            OOFEM_ERROR("Unknown value of damage flag. Must be 0, 1, 2 or 3");
        }
    } else {
        stress = effectiveStress;
    }

    status->letTempStrainVectorBe(fullStrainVector);
    status->letTempAlphaBe(alpha);
    status->letTempStressVectorBe(stress);
    status->letTempEffectiveStressBe(effectiveStress);
#ifdef keep_track_of_dissipated_energy
    double gf = pow(ft, 2) / this->eM; //rough estimation only for this purpose
    status->computeWork(gp, gf);
#endif
    assignStateFlag(gp);
    return stress;
}


FloatArrayF< 2 >
ConcreteDPM2::computeDamage(const FloatArrayF< 6 > &strain,
                            const FloatMatrixF< 6, 6 > &D,
                            double deltaTime,
                            GaussPoint *gp,
                            TimeStep *tStep,
                            double tempAlpha,
                            const FloatArrayF< 6 > &effectiveStress) const
{
    auto status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    double tempEquivStrain;
    double deltaPlasticStrainNorm;
    double tempDamageTension = 0.0;
    double tempDamageCompression = 0.0;

    double tempKappaDTension = 0.0, tempKappaDCompression = 0.0;
    double tempKappaDTensionOne = 0.0, tempKappaDTensionTwo = 0.0;
    double tempKappaDCompressionOne = 0.0, tempKappaDCompressionTwo = 0.0;

    double minEquivStrain = 0.;

    double sig, rho, theta;
    //Calculate coordinates
    computeCoordinates(effectiveStress, sig, rho, theta);

    int unAndReloadingFlag = checkForUnAndReloading(tempEquivStrain, minEquivStrain, D, gp);

    double rateFactor;
    if ( ( status->giveDamageTension() == 0. ) && ( status->giveDamageCompression() == 0. ) ) {
        rateFactor = computeRateFactor(tempAlpha, deltaTime, gp, tStep);
    } else {
        rateFactor = status->giveRateFactor();
    }


    //Compute equivalent strains for  tension and compression
    double tempEquivStrainTension = 0.;
    double tempEquivStrainCompression = 0.;

    tempEquivStrainTension = status->giveEquivStrainTension() + ( tempEquivStrain - status->giveEquivStrain() ) / rateFactor;

    if ( unAndReloadingFlag == 0 ) { //Standard way
        tempEquivStrainCompression = status->giveEquivStrainCompression() + ( tempAlpha * ( tempEquivStrain - status->giveEquivStrain() ) ) / rateFactor;
    } else {
        tempEquivStrainCompression = status->giveEquivStrainCompression() + status->giveAlpha() * ( minEquivStrain - status->giveEquivStrain() ) / rateFactor + ( tempAlpha * ( tempEquivStrain - minEquivStrain ) ) / rateFactor;
    }


    //If damage threshold is exceeded determine the rate factor from the previous step
    if ( ( tempEquivStrainTension > e0 || tempEquivStrainCompression > e0 ) &&
         ( ( status->giveDamageTension() == 0. ) && ( status->giveDamageCompression() == 0. ) ) && !tStep->isTheFirstStep() ) {
        //Rate factor from last step
        rateFactor = status->giveRateFactor();

        tempEquivStrainTension = status->giveEquivStrainTension() + ( tempEquivStrain - status->giveEquivStrain() ) / rateFactor;
        if ( unAndReloadingFlag == 0 ) { //Standard way
            tempEquivStrainCompression = status->giveEquivStrainCompression() + ( tempAlpha * ( tempEquivStrain - status->giveEquivStrain() ) ) / rateFactor;
        } else {
            tempEquivStrainCompression = status->giveEquivStrainCompression() + status->giveAlpha() * ( minEquivStrain - status->giveEquivStrain() ) / rateFactor + ( tempAlpha * ( tempEquivStrain - minEquivStrain ) ) / rateFactor;
        }
    }

    status->letTempRateFactorBe(rateFactor);

    double fTension = tempEquivStrainTension - status->giveKappaDTension();
    double fCompression = tempEquivStrainCompression - status->giveKappaDCompression();

    //Normalize the fs
    fTension = fTension / e0;
    fCompression = fCompression / e0;

    double ductilityMeasure = computeDuctilityMeasureDamage(gp, sig, rho);
    double deltaPlasticStrainNormTension, deltaPlasticStrainNormCompression;

    if ( fTension < -yieldTolDamage && fCompression < -yieldTolDamage ) {
        //Neither tension nor compression is active

        tempKappaDTension = status->giveKappaDTension();
        tempKappaDTensionOne = status->giveKappaDTensionOne();
        tempKappaDTensionTwo = status->giveKappaDTensionTwo();

        tempKappaDCompression = status->giveKappaDCompression();
        tempKappaDCompressionOne = status->giveKappaDCompressionOne();
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo();

        tempDamageTension = status->giveDamageTension();
        tempDamageCompression = status->giveDamageCompression();
    } else if ( fTension >= -yieldTolDamage && fCompression < -yieldTolDamage ) {   //Only tension is active
        //Update tension history variables
        tempKappaDTension = tempEquivStrainTension;
        deltaPlasticStrainNorm = computeDeltaPlasticStrainNormTension(tempKappaDTension, status->giveKappaDTension(), gp);
        tempKappaDTensionOne = status->giveKappaDTensionOne() + deltaPlasticStrainNorm / ductilityMeasure / rateFactor;
        tempKappaDTensionTwo = status->giveKappaDTensionTwo() + ( tempKappaDTension - status->giveKappaDTension() ) / ductilityMeasure;

        //Nothing changes for compression history variables
        tempKappaDCompression = status->giveKappaDCompression();
        tempKappaDCompressionOne = status->giveKappaDCompressionOne();
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo();

        //Initialise damage with tensile history variable
        this->initDamaged(tempKappaDTension, strain, gp);

        tempDamageTension = computeDamageParamTension(tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, status->giveLe(), status->giveDamageTension(), rateFactor);

        tempDamageCompression = status->giveDamageCompression();
    } else if ( fTension < -yieldTolDamage && fCompression >= -yieldTolDamage ) {
        //Only compression is active

        //Nothing changes for the history variables in tension
        tempKappaDTension = status->giveKappaDTension();
        tempKappaDTensionOne = status->giveKappaDTensionOne();
        tempKappaDTensionTwo = status->giveKappaDTensionTwo();

        //Update compression history variables
        tempKappaDCompression = tempEquivStrainCompression;
        deltaPlasticStrainNormCompression = computeDeltaPlasticStrainNormCompression(tempAlpha, tempKappaDCompression, status->giveKappaDCompression(), gp, rho);
        tempKappaDCompressionOne = status->giveKappaDCompressionOne() + deltaPlasticStrainNormCompression / ( ductilityMeasure * rateFactor );
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo() + ( tempKappaDCompression - status->giveKappaDCompression() ) / ductilityMeasure;

        //Determine damage parameters
        tempDamageTension = status->giveDamageTension();
        tempDamageCompression = computeDamageParamCompression(tempKappaDCompression, tempKappaDCompressionOne, tempKappaDCompressionTwo, status->giveDamageCompression(), rateFactor);
    } else if ( fTension >= -yieldTolDamage && fCompression >= -yieldTolDamage ) {
        //Both tension and compression is active

        //Update tension history variables
        tempKappaDTension = tempEquivStrainTension;
        deltaPlasticStrainNormTension = computeDeltaPlasticStrainNormTension(tempKappaDTension, status->giveKappaDTension(), gp);
        tempKappaDTensionOne = status->giveKappaDTensionOne() + deltaPlasticStrainNormTension / ( ductilityMeasure * rateFactor );
        tempKappaDTensionTwo = status->giveKappaDTensionTwo() + ( tempKappaDTension - status->giveKappaDTension() ) / ductilityMeasure;

        //Update the compression history variables
        tempKappaDCompression = tempEquivStrainCompression;
        deltaPlasticStrainNormCompression =
            computeDeltaPlasticStrainNormCompression(tempAlpha, tempKappaDCompression, status->giveKappaDCompression(), gp, rho);
        tempKappaDCompressionOne = status->giveKappaDCompressionOne() + deltaPlasticStrainNormCompression / ( ductilityMeasure * rateFactor );
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo() +
                                   ( tempKappaDCompression - status->giveKappaDCompression() ) / ductilityMeasure;

        //Determine the damage parameters
        this->initDamaged(tempKappaDTension, strain, gp);

        tempDamageTension = computeDamageParamTension(tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, status->giveLe(), status->giveDamageTension(), rateFactor);

        tempDamageCompression = computeDamageParamCompression(tempKappaDCompression, tempKappaDCompressionOne, tempKappaDCompressionTwo, status->giveDamageCompression(), rateFactor);
    }

    //Write all temp history variables to the status
    status->letTempEquivStrainBe(tempEquivStrain);

    //Tension
    status->letTempEquivStrainTensionBe(tempEquivStrainTension);
    status->letTempKappaDTensionBe(tempKappaDTension);
    status->letTempKappaDTensionOneBe(tempKappaDTensionOne);
    status->letTempKappaDTensionTwoBe(tempKappaDTensionTwo);
    status->letTempDamageTensionBe(tempDamageTension);

    //Compression
    status->letTempEquivStrainCompressionBe(tempEquivStrainCompression);
    status->letTempKappaDCompressionBe(tempKappaDCompression);
    status->letTempKappaDCompressionOneBe(tempKappaDCompressionOne);
    status->letTempKappaDCompressionTwoBe(tempKappaDCompressionTwo);
    status->letTempDamageCompressionBe(tempDamageCompression);

    return {
        tempDamageTension, tempDamageCompression
    };
}

int
ConcreteDPM2::checkForUnAndReloading(double &tempEquivStrain,
                                     double &minEquivStrain,
                                     const FloatMatrixF< 6, 6 > &D,
                                     GaussPoint *gp) const
{
    auto status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    //Access old and new strains
    const auto &oldStrain = status->giveReducedStrain();
    const auto &strain = status->giveTempReducedStrain();

    //Compute the temp equivalent strain
    auto tempElasticStrain = strain - status->giveTempPlasticStrain();
    auto tempEffectiveStress = dot(D, tempElasticStrain);

    double sigEffective, rhoEffective, thetaEffective;
    computeCoordinates(tempEffectiveStress, sigEffective, rhoEffective, thetaEffective);
    tempEquivStrain = computeEquivalentStrain(sigEffective, rhoEffective, thetaEffective);
    //Get the equivalent strain from the status
    double equivStrain = status->giveEquivStrain();

    //Compute the increment of effective stress
    auto elasticStrain = oldStrain - status->givePlasticStrain();
    auto effectiveStress = dot(D, elasticStrain);

    auto deltaEffectiveStress = tempEffectiveStress - effectiveStress;

    //Compute equivalent strains for stress state slightly greater than the effective stress and smaller than the temp effective stress
    //For slightly more than effective stress
    auto intermediateEffectiveStressPlus = effectiveStress + 0.01 * deltaEffectiveStress;
    computeCoordinates(intermediateEffectiveStressPlus, sigEffective, rhoEffective, thetaEffective);
    double equivStrainPlus = computeEquivalentStrain(sigEffective, rhoEffective, thetaEffective);

    //For slightly less than temp effective stress
    auto intermediateEffectiveStressMinus = effectiveStress + 0.99 * deltaEffectiveStress;
    computeCoordinates(intermediateEffectiveStressMinus, sigEffective, rhoEffective, thetaEffective);
    double tempEquivStrainMinus = computeEquivalentStrain(sigEffective, rhoEffective, thetaEffective);

    //Check for unloading and reloading in the same step
    int unloadingFlag = 0;
    minEquivStrain = equivStrain;

    if ( ( equivStrain > equivStrainPlus && tempEquivStrain > tempEquivStrainMinus ) &&
         ( fabs(equivStrainPlus - equivStrain) > yieldTolDamage / 100. && fabs(tempEquivStrainMinus - tempEquivStrain) > yieldTolDamage / 100. ) ) {
        unloadingFlag = 1;
        //Unloading and reloading takes place. Find the minimum equivalent strain by subincrementing the effective stress increment
        for ( double k = 1.0; k <= 100.0; k = k + 1.0 ) {
            auto intermediateEffectiveStress = effectiveStress + k / 100. * deltaEffectiveStress;
            computeCoordinates(intermediateEffectiveStress, sigEffective, rhoEffective, thetaEffective);
            double midEquivStrain = computeEquivalentStrain(sigEffective, rhoEffective, thetaEffective);

            if ( midEquivStrain <= minEquivStrain ) {
                minEquivStrain = midEquivStrain;
            } else {
                return unloadingFlag;
            }
        }
    }
    return unloadingFlag;
}


double
ConcreteDPM2::computeRateFactor(double alpha,
                                double deltaTime,
                                GaussPoint *gp,
                                TimeStep *tStep) const
{
    if ( this->strengthRateType == 0 ) {
        return 1;
    }

    auto status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    const auto &strain = status->giveTempReducedStrain();

    //Determine the principal values of the strain
    auto principalStrain = StructuralMaterial::computePrincipalValues(from_voigt_strain(strain) );   ///@todo CHECK

    //Determine max and min value;
    double maxStrain = -1.e20, minStrain = 1.e20;
    for ( int k = 1; k <= principalStrain.giveSize(); k++ ) {
        //maximum
        if ( principalStrain.at(k) > maxStrain ) {
            maxStrain = principalStrain.at(k);
        }

        //minimum
        if ( principalStrain.at(k) < minStrain ) {
            minStrain = principalStrain.at(k);
        }
    }

    //Evaluate the equivalent strains
    double strainRate;
    double oldRateStrain = status->giveRateStrain();
    if ( 1. - alpha > CDPM2_TOL ) { //Tension
        strainRate = ( maxStrain - oldRateStrain ) / deltaTime;
        status->letTempRateStrainBe(maxStrain);
    } else { //Compression
        strainRate = ( minStrain - oldRateStrain ) / deltaTime;
        status->letTempRateStrainBe(minStrain);
    }

    //Tension
    //For tension according to Model Code 2010
    double rateFactorTension = 1.;
    double strainRateRatioTension = strainRate / 1.e-6;

    if ( this->strengthRateType == 1 ) {
        if ( strainRate < 1.e-6 ) {
            rateFactorTension = 1.;
        } else if ( 1.e-6 < strainRate ) {
            rateFactorTension = pow(strainRateRatioTension, 0.018);
        }
    } else if ( this->strengthRateType == 2 ) {
        if ( strainRate < 1.e-6 ) {
            rateFactorTension = 1.;
        } else if ( 1.e-6 < strainRate && strainRate < 10 ) {
            rateFactorTension = pow(strainRateRatioTension, 0.018);
        } else {
            rateFactorTension =  0.0062 * pow(strainRateRatioTension, 1. / 3.);
        }
    }

    //For compression according to Model Code 2010
    double rateFactorCompression = 1.;
    double strainRateRatioCompression = strainRate / ( -30.e-6 );
    if ( this->strengthRateType == 1 ) {
        if ( strainRate > -30.e-6 ) {
            rateFactorCompression = 1.;
        } else if ( -30.e-6 > strainRate ) {
            rateFactorCompression = pow(strainRateRatioCompression, 0.014);
        }
    } else if ( this->strengthRateType == 2 ) {
        if ( strainRate > -30.e-6 ) {
            rateFactorCompression = 1.;
        } else if ( -30.e-6 > strainRate && strainRate > -30 ) {
            rateFactorCompression = pow(strainRateRatioCompression, 0.014);
        } else if ( -30 > strainRate && strengthRateType == 2 ) {
            rateFactorCompression =  0.012 * pow(strainRateRatioCompression, 0.333);
        }
    }

    double rateFactor = ( 1. - alpha ) * rateFactorTension + alpha * rateFactorCompression;

    return rateFactor;
}


double
ConcreteDPM2::computeDeltaPlasticStrainNormTension(double tempKappaD, double kappaD, GaussPoint *gp) const
{
    auto status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    const auto &tempPlasticStrain = status->giveTempPlasticStrain();
    const auto &plasticStrain = status->givePlasticStrain();

    auto deltaPlasticStrain = tempPlasticStrain - plasticStrain;

    double deltaPlasticStrainNorm = 0;

    //Distinguish pre-peak, peak and post-peak

    double factor = 0.;
    if ( tempKappaD < e0 * ( 1. - yieldTolDamage ) ) {
        deltaPlasticStrainNorm = 0.;
    } else if ( tempKappaD > e0 * ( 1. - yieldTolDamage ) && kappaD < e0  * ( 1. - yieldTolDamage ) ) {
        factor = ( 1. - ( e0 - kappaD ) / ( tempKappaD - kappaD ) );
        deltaPlasticStrain *= factor;
        deltaPlasticStrainNorm = norm(deltaPlasticStrain);
    } else {
        deltaPlasticStrainNorm = norm(deltaPlasticStrain);
    }

    return deltaPlasticStrainNorm;
}


double
ConcreteDPM2::computeDeltaPlasticStrainNormCompression(double tempAlpha, double tempKappaD, double kappaD, GaussPoint *gp, const double rho) const
{
    auto status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    const auto &tempPlasticStrain = status->giveTempPlasticStrain();
    const auto &plasticStrain = status->givePlasticStrain();

    auto deltaPlasticStrain = tempAlpha * ( tempPlasticStrain - plasticStrain );

    double deltaPlasticStrainNorm = 0;

    //Distinguish pre-peak, peak and post-peak
    if ( tempKappaD < e0 * ( 1. - yieldTolDamage ) ) {
        deltaPlasticStrainNorm = 0.;
    } else if ( tempKappaD > e0 * ( 1. - yieldTolDamage ) && kappaD < e0 * ( 1. - yieldTolDamage ) ) {
        double factor = ( 1. - ( e0 - kappaD ) / ( tempKappaD - kappaD ) );
        deltaPlasticStrain *= factor;
        deltaPlasticStrainNorm = norm(deltaPlasticStrain);
    } else {
        deltaPlasticStrainNorm = norm(deltaPlasticStrain);
    }

    double tempKappaP = status->giveTempKappaP();
    double yieldHardTwo = computeHardeningTwo(tempKappaP);
    double extraFactor;
    if ( rho < 1.e-16 ) {
        extraFactor = this->ft * yieldHardTwo * sqrt(2. / 3.) / 1.e-16 / sqrt(1. + 2. * pow(this->dilationConst, 2.) );
    } else {
        extraFactor = this->ft * yieldHardTwo * sqrt(2. / 3.) / rho / sqrt(1. + 2. * pow(this->dilationConst, 2.) );
    }

    return deltaPlasticStrainNorm * extraFactor;
}


double
ConcreteDPM2::computeEquivalentStrain(double sig,
                                      double rho,
                                      double theta) const
{
    double rFunction = ( 4. * ( 1. - pow(this->ecc, 2.) ) * pow(cos(theta), 2.) + pow(2. * this->ecc - 1., 2.) ) / ( 2. * ( 1. - pow(this->ecc, 2.) ) * cos(theta) + ( 2. * this->ecc - 1. ) * sqrt(4. * ( 1. - pow(this->ecc, 2.) ) * pow(cos(theta), 2.) + 5. * pow(this->ecc, 2.) - 4. * this->ecc) );

    double pHelp = -this->m * ( rho * rFunction / ( sqrt(6.) * fc ) + sig / fc );

    double qHelp = -3. / 2. * pow(rho, 2.) / pow(this->fc, 2.);

    double help = -0.5 * pHelp + sqrt(pow(pHelp, 2.) / 4. - qHelp);

    double tempEquivStrain = 0.;
    if ( help > 0 ) {
        tempEquivStrain = help * e0;
    }
    return tempEquivStrain;
}


double
ConcreteDPM2::computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld, double rateFactor) const
{
    double omega = 0.;

    //So that damage does not turn out to be negative if function is entered for equivstrains smaller thatn e0.
    double ftTemp = this->ft * ( 1. - yieldTolDamage );

    double wfMod = this->wf;
    double wfOneMod = this->wfOne;

    if ( this->strengthRateType > 0 ) {
        if ( this->energyRateType == 0 ) {
            wfMod /= pow(rateFactor, 2.);
            wfOneMod /= pow(rateFactor, 2.);
        } else if ( this->energyRateType == 1 ) {
            wfMod /= rateFactor;
            wfOneMod /= rateFactor;
        }
    }

    double help;
    if ( equivStrain > e0 * ( 1. - yieldTolDamage ) ) {
        if ( softeningType == 0 ) { //linear
            omega = ( this->eM * equivStrain * wfMod - ftTemp * wfMod + ftTemp * kappaOne * le ) /
                    ( this->eM * equivStrain * wfMod - ftTemp * le * kappaTwo );
        } else if ( softeningType == 1 ) { //bilinear: Calculate damage parameter for both parts of bilinear curve  and check which fulfils limits.
            omega = ( this->eM * equivStrain * wfOneMod - ftTemp * wfOneMod - ( this->ftOne - ftTemp ) * kappaOne * le ) /
                    ( this->eM * equivStrain * wfOneMod + ( this->ftOne - ftTemp ) * le * kappaTwo );
            help = le * kappaOne + le * omega * kappaTwo;

            if ( help >= 0. && help < wfOneMod ) {
                return omega;
            }

            omega = ( this->eM * equivStrain * ( wfMod - wfOneMod ) - this->ftOne * ( wfMod - wfOneMod ) +
                      this->ftOne * kappaOne * le  - this->ftOne * wfOneMod ) /
                    ( this->eM * equivStrain * ( wfMod - wfOneMod )  - this->ftOne * le * kappaTwo );
            help = le * kappaOne + le * omega * kappaTwo;

            if ( help > wfOneMod && help < wfMod ) {
                return omega;
            }
        } else if ( softeningType == 2 ) { //exponential: Iterative solution
            omega = 1.; //Initial guess
            double residual = 0.;
            double dResidualDOmega = 0.;
            int nite = 0;

            do {
                nite++;

                residual  = ( 1 - omega ) * this->eM * equivStrain - ftTemp * exp(-le * ( omega * kappaTwo + kappaOne ) / wfMod);
                dResidualDOmega = -this->eM * equivStrain + ftTemp * le * kappaTwo / wfMod * exp(-le * ( omega * kappaTwo + kappaOne ) / wfMod);

                omega -= residual / dResidualDOmega;
                if ( nite > newtonIter ) {
                    OOFEM_ERROR("algorithm not converging");
                }
            } while ( fabs(residual / ftTemp) >= 1.e-8 );
        }
    } else {
        omega = 0.;
    }


    if ( omega > 1. ) {
        omega = 1.;
    }

    if ( omega < 0. || omega < omegaOld ) {
        omega = omegaOld;
    }


    return omega;
}

double
ConcreteDPM2::computeDamageParamCompression(double equivStrain, double kappaOne, double kappaTwo, double omegaOld, double rateFactor) const
{
    if ( this->damageFlag == 3 ) {
        return 0.;
    }

    double ftTemp = this->ft * ( 1. - yieldTolDamage );
    double efCompressionMod = this->efCompression;

    if ( this->strengthRateType > 0 ) {
        if ( this->energyRateType == 0 ) {
            efCompressionMod /= pow(rateFactor, 2.);
        } else if ( this->energyRateType == 1 ) {
            efCompressionMod /= rateFactor;
        }
    }

    double omega = 1.;
    int nite = 0;
    double residual = 0.;
    double dResidualDOmega = 0.;

    if ( equivStrain > e0 * ( 1. - yieldTolDamage ) ) {
        do {
            nite++;

            residual = ( 1. - omega ) * this->eM * equivStrain - ftTemp * exp(-( kappaOne + omega * kappaTwo ) / efCompressionMod);
            dResidualDOmega = -this->eM * equivStrain + ftTemp * kappaTwo / efCompressionMod * exp(-( kappaOne + omega * kappaTwo ) / efCompressionMod);

            omega -= residual / dResidualDOmega;
            if ( nite > newtonIter ) {
                OOFEM_ERROR("algorithm not converging");
            }
        } while ( fabs(residual / ft) >= 1.e-8 );
    } else {
        omega = 0.;
    }

    if ( omega > 1. ) {
        omega = 1.;
    }
    if ( omega < omegaOld || omega < 0. ) {
        omega = omegaOld;
    }

    return omega;
}



void
ConcreteDPM2::initDamaged(double kappaD,
                          const FloatArrayF< 6 > &strain,
                          GaussPoint *gp) const
{
    if ( kappaD <= e0 * ( 1. - yieldTolDamage ) ) {
        return;
    }

    auto status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    if ( helem > 0. ) {
        status->setLe(helem);
    } else if ( status->giveDamageTension() == 0. && status->giveDamageCompression() == 0. ) {
        //auto [principalStrains, principalDir] = computePrincipalValDir(from_voigt_strain(strain)); // c++17
        auto tmp = computePrincipalValDir(from_voigt_strain(strain) );
        auto principalStrains = tmp.first;
        auto principalDir = tmp.second;

        // find index of max positive principal strain
        int indx = 1;
        for ( int i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) > principalStrains.at(indx) ) {
                indx = i;
            }
        }

        FloatArray crackPlaneNormal(3);
        for ( int i = 1; i <= 3; i++ ) {
            crackPlaneNormal.at(i) = principalDir.at(i, indx);
        }

        // evaluate the projected element size
        double le = gp->giveElement()->giveCharacteristicLength(crackPlaneNormal);
        if ( le == 0. ) {
            le = gp->giveElement()->computeMeanSize();
        }

        // store le in the corresponding status
        status->setLe(le);
    } else if ( status->giveLe() == 0. ) {
        // this happens if the status is initialized from a file
        // with nonzero damage
        // le determined as square root of element area or cube root of el. volume
        double le = gp->giveElement()->computeMeanSize();
        status->setLe(le);
    }
}


double
ConcreteDPM2::computeDuctilityMeasureDamage(GaussPoint *gp, const double sig, const double rho) const
{
    //1./sqrt(6.)=0.40824829
    double alphaZero = 0.40824829;

    double Rs = 0;
    if ( sig < 0. ) {
        if ( rho > 1.e-16 ) {
            Rs = -sig / ( alphaZero * rho );
        } else { //Need to set a dummy valye
            Rs = -sig * 1.e16 / alphaZero;
        }
    } else {
        Rs = 0;
    }

    return 1. + ( ASoft - 1. ) * Rs; // ductilityMeasure
}


FloatArrayF< 6 >
ConcreteDPM2::performPlasticityReturn(GaussPoint *gp, const FloatMatrixF< 6, 6 > &D, const FloatArrayF< 6 > &strain) const
{
    auto status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    //get plastic strain and kappa
    auto tempPlasticStrain = status->givePlasticStrain();
    double tempKappaP = status->giveKappaP();

    //this theta computed here should stay constant for the rest of procedure.

    const auto &oldStrain = status->giveReducedStrain();

    // introduce a strange subincrementation flag
    int subIncrementFlag = 0;

    double apexStress = 0.;
    int subincrementcounter = 0;
    //Attempt to implement subincrementation
    // initialize variables
    subIncrementFlag = 0;
    auto convergedStrain = oldStrain;
    auto tempStrain = strain;
    auto deltaStrain = strain - oldStrain;

    FloatArrayF< 6 >effectiveStress;

    //To get into the loop
    status->letTempReturnResultBe(ConcreteDPM2Status::RR_NotConverged);
    while ( status->giveTempReturnResult() == ConcreteDPM2Status::RR_NotConverged || subIncrementFlag == 1 ) {
        auto elasticStrain = tempStrain - tempPlasticStrain;

        effectiveStress = dot(D, elasticStrain);

        double sig, rho, theta;
        computeCoordinates(effectiveStress, sig, rho, theta);
        double yieldValue = computeYieldValue(sig, rho, theta, tempKappaP);

        apexStress = 0.;

        if ( yieldValue > 0. ) {
            checkForVertexCase(apexStress, sig, tempKappaP, gp);
            if ( status->giveTempReturnType() == ConcreteDPM2Status::RT_Tension || status->giveTempReturnType() == ConcreteDPM2Status::RT_Compression ) {
                tempKappaP = performVertexReturn(effectiveStress, apexStress, tempKappaP, gp);
                status->letTempKappaPBe(tempKappaP);
            }
            if ( status->giveTempReturnType() == ConcreteDPM2Status::RT_Regular ) {
                tempKappaP = performRegularReturn(effectiveStress, tempKappaP, gp, theta);
                status->letTempKappaPBe(tempKappaP);
            }
        } else {
            status->letTempReturnResultBe(ConcreteDPM2Status::RR_Converged);
            tempPlasticStrain = status->givePlasticStrain();
            status->letTempPlasticStrainBe(tempPlasticStrain);
            status->letTempKappaPBe(tempKappaP);
            break;
        }

        if ( status->giveTempReturnResult() == ConcreteDPM2Status::RR_NotConverged ) {
            subincrementcounter++;
            if ( subincrementcounter > 10 ) {
                OOFEM_LOG_INFO("Unstable element %d \n", gp->giveElement()->giveGlobalNumber() );
                OOFEM_LOG_INFO("Old strain vector %g %g %g %g %g %g  \n", oldStrain.at(1), oldStrain.at(2), oldStrain.at(3), oldStrain.at(4), oldStrain.at(5), oldStrain.at(6) );

                const auto &help = status->giveTempPlasticStrain();
                OOFEM_LOG_INFO("Old plastic strain vector %g %g %g %g %g %g  \n", help.at(1), help.at(2), help.at(3), help.at(4), help.at(5), help.at(6) );
                OOFEM_LOG_INFO("New strain vector %g %g %g %g %g %g  \n", strain.at(1), strain.at(2), strain.at(3), strain.at(4), strain.at(5), strain.at(6) );

                computeCoordinates(effectiveStress, sig, rho, theta);
                double sig1, rho1, theta1;
                auto help1 = dot(D, oldStrain - help);
                computeCoordinates(help1, sig1, rho1, theta1);
                yieldValue = computeYieldValue(sig, rho, theta, tempKappaP);
                OOFEM_LOG_INFO("OLD Sig %g rho %g theta %g  \n", sig1, rho1, theta1);
                OOFEM_LOG_INFO("NEW Sig %g rho %g theta %g  \n", sig, rho, theta);
                if ( status->giveTempReturnType() == ConcreteDPM2Status::RT_Tension || status->giveTempReturnType() == ConcreteDPM2Status::RT_Compression ) {
                    OOFEM_LOG_INFO("Vertex case apexstress %g\n", apexStress);
                } else {
                    OOFEM_LOG_INFO("Regular case %g \n", 15.18);
                }
                OOFEM_LOG_INFO("KappaP old %g new %g yieldfun %g\n", status->giveTempKappaP(), tempKappaP, yieldValue);
                OOFEM_ERROR("Could not reach convergence with small deltaStrain, giving up.");
            } else if ( subincrementcounter > 9 && tempKappaP < 1. ) {
                tempKappaP = 1.;
                status->letTempKappaPBe(tempKappaP);
            }

            subIncrementFlag = 1;
            deltaStrain *= 0.5;
            tempStrain = convergedStrain + deltaStrain;
        } else if ( status->giveTempReturnResult() == ConcreteDPM2Status::RR_Converged && subIncrementFlag == 0 ) {
            auto C = inv(D); // compliance
            elasticStrain = dot(C, effectiveStress);
            tempPlasticStrain = strain - elasticStrain;
            status->letTempPlasticStrainBe(tempPlasticStrain);
        } else if ( status->giveTempReturnResult() == ConcreteDPM2Status::RR_Converged && subIncrementFlag == 1 ) {
            subincrementcounter = 0;
            auto C = inv(D); // compliance
            elasticStrain = dot(C, effectiveStress);
            tempPlasticStrain = tempStrain - elasticStrain;
            status->letTempPlasticStrainBe(tempPlasticStrain);

            subIncrementFlag = 0;
            status->letTempReturnResultBe(ConcreteDPM2Status::RR_NotConverged);
            convergedStrain = tempStrain;
            deltaStrain = strain - convergedStrain;
            tempStrain = strain;
        }
    }

    return effectiveStress;
}

void
ConcreteDPM2::checkForVertexCase(double &answer,
                                 double sig,
                                 double tempKappa,
                                 GaussPoint *gp) const
{
    auto status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    answer = 0.;
    if ( sig > 0. ) {
        status->letTempReturnTypeBe(ConcreteDPM2Status::RT_Tension);
    } else if ( sig < 0. &&  tempKappa < 1. ) {
        status->letTempReturnTypeBe(ConcreteDPM2Status::RT_Compression);
    } else {
        status->letTempReturnTypeBe(ConcreteDPM2Status::RT_Regular);
    }
}


double
ConcreteDPM2::performVertexReturn(FloatArrayF< 6 > &effectiveStress,
                                  double apexStress, double tempKappaP,
                                  GaussPoint *gp) const
{
    //auto [deviatoricStressTrial, sigTrial] = computeDeviatoricVolumetricSplit(effectiveStress); // c++17
    auto tmp = computeDeviatoricVolumetricSplit(effectiveStress);
    auto deviatoricStressTrial = tmp.first;
    auto sigTrial = tmp.second;

    auto status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    double rhoTrial = computeSecondCoordinate(deviatoricStressTrial);

    double kappaInitial = tempKappaP;

    double sig2 = apexStress;

    tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sigTrial);

    double yieldValue = computeYieldValue(sigTrial, 0., 0., tempKappaP);

    tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sig2);

    double yieldValueMid = computeYieldValue(sig2, 0., 0., tempKappaP);

    if ( yieldValue * yieldValueMid >= 0. ) {
        status->letTempReturnTypeBe(ConcreteDPM2Status::RT_Regular);
        status->letTempReturnResultBe(ConcreteDPM2Status::RR_NotConverged);
        return kappaInitial;
    }

    double dSig, sigAnswer;
    if ( yieldValue < 0.0 ) {
        dSig = sig2 - sigTrial;
        sigAnswer = sig2;
    } else {
        dSig = sigTrial - sig2;
        sigAnswer = sig2;
    }

    for ( int j = 0; j < 250; j++ ) {
        dSig = 0.5 * dSig;

        double sigMid = sigAnswer + dSig;


        tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sigMid);

        yieldValueMid = computeYieldValue(sigMid, 0., 0., tempKappaP);

        if ( yieldValueMid <= 0. ) {
            sigAnswer = sigMid;
        }

        if ( fabs(yieldValueMid) < yieldTol && yieldValueMid <= 0. ) {
            double ratioPotential =
                computeRatioPotential(sigAnswer,  0, tempKappaP);


            double ratioTrial = rhoTrial / ( sigTrial - sigAnswer );

            if ( ( ( ( ratioPotential >= ratioTrial ) && status->giveTempReturnType() == ConcreteDPM2Status::RT_Tension ) ) ||
                 ( ( ratioPotential <= ratioTrial ) && status->giveTempReturnType() == ConcreteDPM2Status::RT_Compression ) ) {
                for ( int i = 0; i < 3; i++ ) {
                    effectiveStress.at(i + 1) = sigAnswer;
                }

                for ( int i = 3; i < 6; i++ ) {
                    effectiveStress.at(i + 1) = 0.;
                }
                status->letTempReturnResultBe(ConcreteDPM2Status::RR_Converged);
                return tempKappaP;
            } else {
                status->letTempReturnTypeBe(ConcreteDPM2Status::RT_Regular);
                status->letTempReturnResultBe(ConcreteDPM2Status::RR_NotConverged);
                return kappaInitial;
            }
        }
    }

    for ( int i = 0; i < 3; i++ ) {
        effectiveStress.at(i + 1) = sigAnswer;
    }

    for ( int i = 3; i < 6; i++ ) {
        effectiveStress.at(i + 1) = 0.;
    }
    status->letTempReturnResultBe(ConcreteDPM2Status::RR_Converged);

    OOFEM_WARNING("Perform vertex return not converged!\n");

    return tempKappaP;
}


double
ConcreteDPM2::computeTempKappa(double kappaInitial,
                               double sigTrial,
                               double rhoTrial,
                               double sig) const
{
    //This function is called, if stress state is in vertex case
    double equivalentDeltaPlasticStrain = sqrt(1. / 9. * pow( ( sigTrial - sig ) / ( kM ), 2. ) +
                                               pow(rhoTrial / ( 2. * gM ), 2.) );

    double thetaVertex = M_PI / 3.;
    double ductilityMeasure = computeDuctilityMeasure(sig, 0., thetaVertex);

    return kappaInitial + equivalentDeltaPlasticStrain / ductilityMeasure;
}


double
ConcreteDPM2::computeDuctilityMeasure(double sig,
                                      double rho,
                                      double theta) const
{
    double thetaConst = pow(2. * cos(theta), 2.);
    double ductilityMeasure;
    double x = -( sig + fc / 3 ) / fc;
    if ( x < 0. ) {
        /*Introduce exponential help function which results in a smooth
         * transition. */
        double EHard = BHard - DHard;
        double FHard = ( BHard - DHard ) * CHard / ( AHard - BHard );
        ductilityMeasure = ( EHard * exp(x / FHard) + DHard ) / thetaConst;
    } else {
        ductilityMeasure = ( AHard + ( BHard - AHard ) * exp(-x / ( CHard ) ) ) / thetaConst;
    }

    return ductilityMeasure;
}


double
ConcreteDPM2::computeRatioPotential(double sig,
                                    double rho,
                                    double tempKappa) const
{
    //compute yieldHard and yieldSoft
    double yieldHardOne = computeHardeningOne(tempKappa);
    double yieldHardTwo = computeHardeningTwo(tempKappa);

    //Compute dilation parameter
    double AGParam = this->ft * yieldHardTwo * 3 / this->fc + m / 2;
    double BGParam =
        yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) /
        ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

    double R = ( sig - ft / 3. * yieldHardTwo ) / fc / BGParam;
    double mQ = AGParam * exp(R);

    double Bl = sig / fc + rho / ( fc * sqrt(6.) );

    double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fc;

    double dgdsig = 4. * ( 1. - yieldHardOne ) / fc * Al * Bl + pow(yieldHardOne, 2.) * mQ / fc;

    double dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +

                    m * pow(yieldHardOne, 2.) / ( sqrt(6.) * fc );

    return dgdrho / dgdsig * 3. * ( 1. - 2. * nu ) / ( 1. + nu );
}


double
ConcreteDPM2::performRegularReturn(FloatArrayF< 6 > &effectiveStress,
                                   double kappaP,
                                   GaussPoint *gp,
                                   double theta) const
{
    auto status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    //Define stressVariables
    double trialSig, trialRho;

    auto trialStress = effectiveStress;

    //compute invariants from stress state
    //auto [deviatoricTrialStress, trialSig] = computeDeviatoricVolumetricSplit(trialStress); // c++17
    auto tmp = computeDeviatoricVolumetricSplit(trialStress);
    auto deviatoricTrialStress = tmp.first;
    trialSig = tmp.second;
    trialRho = computeSecondCoordinate(deviatoricTrialStress);

    double sig = trialSig;
    double rho = trialRho;

    // Starting guess:
    double tempKappaP = kappaP;

    // Look at the magnitudes of the residuals. You have to scale the yieldValue down.
    double yieldValue = computeYieldValue(sig, rho, theta, tempKappaP);

    //initialise unknowns
    FloatArray unknowns;
    FloatArray residuals;
    residuals.resize(4);
    residuals.at(4) = yieldValue;  //store in the last element of the array
    unknowns.resize(4);
    unknowns.at(1) = trialSig;
    unknowns.at(2) = trialRho;
    unknowns.at(3) = tempKappaP;
    unknowns.at(4) = 0.;

    double deltaLambda = 0.;
    double normOfResiduals  = 1.; //just to get into the loop

    // N.R. iteration for finding the correct plastic return which is found when the norm of the residuals are equal to zero

    int iterationCount = 0;
    while ( normOfResiduals > yieldTol   ) {
        iterationCount++;
        if ( iterationCount == newtonIter ) {
            status->letTempReturnResultBe(ConcreteDPM2Status::RR_NotConverged);

            return kappaP;
        }

        auto residualsNorm = residuals;
        //Normalize residuals. Think about it more.
        residualsNorm.at(1) /= this->kM;
        residualsNorm.at(2) /= 2. * this->gM;

        normOfResiduals = norm(residualsNorm);

        if ( std::isnan(normOfResiduals) ) {
            status->letTempReturnResultBe(ConcreteDPM2Status::RR_NotConverged);
            return kappaP;
        }

        if ( normOfResiduals > yieldTol  ) {
            // Test to run newton iteration using inverse of Jacobian
            auto jacobian = computeJacobian(sig, rho, theta, tempKappaP, deltaLambda, gp);

            try {
                auto deltaIncrement = solve( jacobian, FloatArrayF< 4 >(residuals) );
                unknowns -= deltaIncrement;
            } catch ( ... ) {
                status->letTempReturnResultBe(ConcreteDPM2Status::RR_NotConverged);
                return kappaP;
            }

            unknowns.at(2) = max(unknowns.at(2), 0.); //Keep rho greater than zero!
            unknowns.at(3) = max(unknowns.at(3), kappaP); //Keep deltaKappa greater than zero!
            unknowns.at(4) = max(unknowns.at(4), 0.); //Keep deltaLambda greater than zero!

            //compute residuals
            sig = unknowns.at(1);
            rho = unknowns.at(2);
            tempKappaP = unknowns.at(3);
            deltaLambda = unknowns.at(4);

            /* Compute the mVector holding the derivatives of the g function and the hardening function*/
            auto dGDInv = computeDGDInv(sig, rho, tempKappaP);
            double dKappaDDeltaLambda = computeDKappaDDeltaLambda(sig, rho, theta, tempKappaP);

            residuals.at(1) = sig - trialSig + this->kM * deltaLambda * dGDInv.at(1);
            residuals.at(2) = rho - trialRho + ( 2. * this->gM ) * deltaLambda * dGDInv.at(2);
            residuals.at(3) = -tempKappaP + kappaP + deltaLambda * dKappaDDeltaLambda;
            residuals.at(4) = computeYieldValue(sig, rho, theta, tempKappaP);
        }
    }


    //compute the principal directions of the stress
    //auto [helpStress, stressPrincipalDir] = StructuralMaterial :: computePrincipalValDir(from_voigt_stress(trialStress)); // c++17
    auto tmpEig = StructuralMaterial::computePrincipalValDir( from_voigt_stress(trialStress) );
    auto stressPrincipalDir = tmpEig.second;

    FloatArrayF< 6 >stressPrincipal;
    stressPrincipal [ 0 ] = sig + sqrt(2. / 3.) * rho * cos(theta);
    stressPrincipal [ 1 ] = sig + sqrt(2. / 3.) * rho * cos(theta - 2. * M_PI / 3.);
    stressPrincipal [ 2 ] = sig + sqrt(2. / 3.) * rho * cos(theta + 2. * M_PI / 3.);
    effectiveStress = transformStressVectorTo(stressPrincipalDir, stressPrincipal, 1);
    status->letTempReturnResultBe(ConcreteDPM2Status::RR_Converged);

    //Store deltaLambda in status
    status->letDeltaLambdaBe(deltaLambda);

    return tempKappaP;
}


FloatMatrixF< 4, 4 >
ConcreteDPM2::computeJacobian(double sig,
                              double rho,
                              double theta,
                              double kappa,
                              double deltaLambda,
                              GaussPoint *gp) const
{
    auto dFDInv = computeDFDInv(sig, rho, theta, kappa);
    auto dGDInv = computeDGDInv(sig, rho, kappa);
    auto dDGDDInv = computeDDGDDInv(sig, rho, kappa);

    double dKappaDDeltaLambda = computeDKappaDDeltaLambda(sig, rho, theta, kappa);
    double dFDKappa = computeDFDKappa(sig, rho, theta, kappa);

    auto dDGDInvDKappa = computeDDGDInvDKappa(sig, rho, kappa);

    double dDKappaDDeltaLambdaDKappa = computeDDKappaDDeltaLambdaDKappa(sig, rho, theta, kappa);
    auto dDKappaDDeltaLambdaDInv = computeDDKappaDDeltaLambdaDInv(sig, rho, theta, kappa);

    FloatMatrixF< 4, 4 >answer;
    /* Compute matrix*/
    answer.at(1, 1) = 1. + this->kM * deltaLambda * dDGDDInv.at(1, 1);
    answer.at(1, 2) = this->kM * deltaLambda * dDGDDInv.at(1, 2);
    answer.at(1, 3) = this->kM * deltaLambda * dDGDInvDKappa.at(1);
    answer.at(1, 4) = this->kM * dGDInv.at(1);
    /**/
    answer.at(2, 1) = 2. * this->gM * deltaLambda * dDGDDInv.at(2, 1);
    answer.at(2, 2) = 1. + 2. * this->gM * deltaLambda * dDGDDInv.at(2, 2);
    answer.at(2, 3) = 2. * this->gM * deltaLambda * dDGDInvDKappa.at(2);
    answer.at(2, 4) = 2. * this->gM * dGDInv.at(2);
    /**/
    answer.at(3, 1) = deltaLambda * dDKappaDDeltaLambdaDInv.at(1);
    answer.at(3, 2) = deltaLambda * dDKappaDDeltaLambdaDInv.at(2);
    answer.at(3, 3) = deltaLambda * dDKappaDDeltaLambdaDKappa - 1.;
    answer.at(3, 4) = dKappaDDeltaLambda;
    /**/
    answer.at(4, 1) = dFDInv.at(1);
    answer.at(4, 2) = dFDInv.at(2);
    answer.at(4, 3) = dFDKappa;
    answer.at(4, 4) = 0.;
    return answer;
}




double
ConcreteDPM2::computeYieldValue(double sig,
                                double rho,
                                double theta,
                                double tempKappa) const
{
    //compute yieldHard
    double yieldHardOne = computeHardeningOne(tempKappa);
    double yieldHardTwo = computeHardeningTwo(tempKappa);


    //  compute elliptic function r
    double rFunction = ( 4. * ( 1. - pow(ecc, 2.) ) * pow(cos(theta), 2.) +
                         pow( ( 2. * ecc - 1. ), 2. ) ) /
                       ( 2. * ( 1. - pow(ecc, 2.) ) * cos(theta) +
                         ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - pow(ecc, 2.) ) * pow(cos(theta), 2.)
                                                  + 5. * pow(ecc, 2.) - 4. * ecc) );

    //compute help function Al
    double Al = ( 1. - yieldHardOne ) * pow( ( sig / fc + rho / ( sqrt(6.) * fc ) ), 2. ) +
                sqrt(3. / 2.) * rho / fc;

    //Compute yield equation
    return pow(Al, 2.) +
           pow(yieldHardOne, 2.) * yieldHardTwo * m * ( sig / fc + rho * rFunction / ( sqrt(6.) * fc ) ) -
           pow(yieldHardOne, 2.) * pow(yieldHardTwo, 2.);
}


double
ConcreteDPM2::computeDFDKappa(double sig,
                              double rho,
                              double theta,
                              double tempKappa) const
{
    double dFDKappa;
    //compute yieldHard and yieldSoft
    double yieldHardOne = computeHardeningOne(tempKappa);
    double yieldHardTwo = computeHardeningTwo(tempKappa);
    // compute the derivative of the hardening and softening laws
    double dYieldHardOneDKappa = computeHardeningOnePrime(tempKappa);
    double dYieldHardTwoDKappa = computeHardeningTwoPrime(tempKappa);
    //compute elliptic function r
    double rFunction =
        ( 4. * ( 1. - pow(ecc, 2) ) * pow(cos(theta), 2) + pow( ( 2. * ecc - 1. ), 2 ) ) /
        ( 2 * ( 1. - pow(ecc, 2) ) * cos(theta) + ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - pow(ecc, 2) ) * pow(cos(theta), 2) + 5. * pow(ecc, 2) - 4. * ecc) );

    //compute help functions Al, Bl
    double Al = ( 1. - yieldHardOne ) * pow( ( sig / fc + rho / ( sqrt(6.) * fc ) ), 2.) + sqrt(3. / 2.) * rho / fc;


    double Bl = sig / fc + rho / ( fc * sqrt(6.) );
    double dFDYieldHardOne = -2. * Al * pow(Bl, 2.)
                             + 2. * yieldHardOne * yieldHardTwo * m * ( sig / fc + rho * rFunction / ( sqrt(6.) * fc ) ) - 2. * yieldHardOne * pow(yieldHardTwo, 2.);

    double dFDYieldHardTwo = pow(yieldHardOne, 2.) * m * ( sig / fc + rho * rFunction / ( sqrt(6.) * fc ) ) - 2. * yieldHardTwo * pow(yieldHardOne, 2.);

    // compute dFDKappa
    dFDKappa = dFDYieldHardOne * dYieldHardOneDKappa +
               dFDYieldHardTwo * dYieldHardTwoDKappa;
    /*
     * set dFDKappa to zero, if it becomes greater than zero.
     * dFDKappa can only be negative or zero in the converged state for
     * the case of hardenig and perfect plasticity. For trial stresses, however,
     * it might be psoitive, which may lead to convergency problems. Therefore,
     * it is set to zero in this cases.
     */
    if ( dFDKappa > 0. ) {
        dFDKappa = 0.;
    }

    return dFDKappa;
}

FloatArrayF< 2 >
ConcreteDPM2::computeDFDInv(double sig,
                            double rho,
                            double theta,
                            double tempKappa) const
{
    //compute yieldHard
    double yieldHardOne = computeHardeningOne(tempKappa);
    double yieldHardTwo = computeHardeningTwo(tempKappa);

    //compute elliptic function r
    double rFunction = ( 4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta) + ( 2. * ecc - 1. ) * ( 2. * ecc - 1. ) ) /
                       ( 2. * ( 1. - ecc * ecc ) * cos(theta) + ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta)
                                                                                         + 5. * ecc * ecc - 4. * ecc) );

    //compute help functions AL, BL
    double AL = ( 1. - yieldHardOne ) * pow( ( sig / fc + rho / ( sqrt(6.) * fc ) ), 2. ) + sqrt(3. / 2.) * rho / fc;
    double BL = sig / fc + rho / ( fc * sqrt(6.) );

    //compute dfdsig
    double dfdsig = 4. * ( 1. - yieldHardOne ) / fc * AL * BL + yieldHardTwo * pow(yieldHardOne, 2.) * m / fc;

    //compute dfdrho
    double dfdrho = AL / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * BL + 6. ) + rFunction * m * yieldHardTwo * pow(yieldHardOne, 2.) / ( sqrt(6.) * fc );

    return {
        dfdsig, dfdrho
    };
}


double
ConcreteDPM2::computeDKappaDDeltaLambda(double sig,
                                        double rho,
                                        double theta,
                                        double tempKappa) const
{
    auto dGDInv = computeDGDInv(sig, rho, tempKappa);
    double equivalentDGDStress = sqrt(1. / 3. * pow(dGDInv [ 0 ], 2.) + pow(dGDInv [ 1 ], 2.) );
    double ductilityMeasure = computeDuctilityMeasure(sig, rho, theta);
    return equivalentDGDStress / ductilityMeasure; // dKappaDDeltaLambda
}

FloatArrayF< 2 >
ConcreteDPM2::computeDDKappaDDeltaLambdaDInv(double sig,
                                             double rho,
                                             double theta,
                                             double tempKappa) const
{
    //Compute first and second derivative of plastic potential
    auto dGDInv = computeDGDInv(sig, rho, tempKappa);
    auto dDGDDInv = computeDDGDDInv(sig, rho, tempKappa);

    //Compute equivalentDGDStress
    double equivalentDGDStress = sqrt(1. / 3. * pow(dGDInv [ 0 ], 2.) + pow(dGDInv [ 1 ], 2.) );

    //computeDuctilityMeasure
    double ductilityMeasure = computeDuctilityMeasure(sig, rho, theta);

    //Compute dEquivalentDGDStressDInv
    FloatArrayF< 2 >dEquivalentDGDStressDInv;
    dEquivalentDGDStressDInv [ 0 ] =
        ( 2. / 3. * dGDInv [ 0 ] * dDGDDInv(0, 0) + 2. * dGDInv [ 1 ] * dDGDDInv(1, 0) ) / ( 2. * equivalentDGDStress );
    dEquivalentDGDStressDInv [ 1 ] =
        ( 2. / 3. * dGDInv [ 0 ] * dDGDDInv(0, 1) + 2. * dGDInv [ 1 ] * dDGDDInv(1, 1) ) / ( 2. * equivalentDGDStress );


    // compute the derivative of
    auto dDuctilityMeasureDInv = computeDDuctilityMeasureDInv(sig, rho, theta, tempKappa);

    FloatArrayF< 2 >answer;
    answer [ 0 ] = ( dEquivalentDGDStressDInv [ 0 ] * ductilityMeasure - equivalentDGDStress * dDuctilityMeasureDInv [ 0 ] ) / pow(ductilityMeasure, 2.);
    answer [ 1 ] = ( dEquivalentDGDStressDInv [ 1 ] * ductilityMeasure - equivalentDGDStress * dDuctilityMeasureDInv [ 1 ] ) / pow(ductilityMeasure, 2.);
    return answer;
}

FloatArrayF< 6 >
ConcreteDPM2::computeDDKappaDDeltaLambdaDStress(const FloatArrayF< 6 > &stress, double tempKappa) const
{
    auto tmp = computeDeviatoricVolumetricSplit(stress);
    auto deviatoricStress = tmp.first;
    double sig = tmp.second;

    double rho = computeSecondCoordinate(deviatoricStress);
    double theta = computeThirdCoordinate(deviatoricStress);

    auto dDKappaDDeltaLambdaDInv = computeDDKappaDDeltaLambdaDInv(sig, rho, theta, tempKappa);

    // compute dDKappaDDeltaLambdaDCosTheta
    auto dGDInv = computeDGDInv(sig, rho, tempKappa);

    double equivalentDGDStress = sqrt(1. / 3. * pow(dGDInv [ 0 ], 2.) + pow(dGDInv [ 1 ], 2.) );

    double ductilityMeasure = computeDuctilityMeasure(sig, rho, theta);

    //Reuse implementation to compute dKappaDDeltaLambdaDCosTheta
    //Ductility measure has in denominator (2*cos(theta))^2
    double dDKappaDDeltaLambdaDCosTheta = equivalentDGDStress / ductilityMeasure / cos(theta);
    //  dDKappaDDeltaLambdaDCosTheta = 0.;

    auto dSigDStress = computeDSigDStress();
    auto dRhoDStress = computeDRhoDStress(stress);
    auto dCosThetaDStress = computeDCosThetaDStress(stress);

    dSigDStress *= dDKappaDDeltaLambdaDInv.at(1);

    dRhoDStress *= dDKappaDDeltaLambdaDInv.at(2);

    dCosThetaDStress *= dDKappaDDeltaLambdaDCosTheta;

    auto dDKappaDDeltaLambdaDStress = dSigDStress + dRhoDStress + dCosThetaDStress;

    return dDKappaDDeltaLambdaDStress;
}

FloatArrayF< 6 >
ConcreteDPM2::computeDDGDStressDKappa(const FloatArrayF< 6 > &stress, double tempKappa) const
{
    FloatArrayF< 6 >answer;

    auto tmp = computeDeviatoricVolumetricSplit(stress);
    auto deviatoricStress = tmp.first;
    double sig = tmp.second;

    double rho = computeSecondCoordinate(deviatoricStress);
    //  double theta = computeThirdCoordinate(deviatoricStress);

    auto dDGDInvDKappa = computeDDGDInvDKappa(sig, rho, tempKappa);
    auto dSigDStress = computeDSigDStress();
    auto dRhoDStress = computeDRhoDStress(stress);

    answer = dSigDStress;
    answer *= dDGDInvDKappa.at(1);

    FloatArrayF< 6 >temp1;
    temp1 = dRhoDStress;
    temp1 *= dDGDInvDKappa.at(2);

    answer += temp1;

    return answer;
}



double
ConcreteDPM2::computeDDKappaDDeltaLambdaDKappa(double sig, double rho, double theta, double tempKappa) const
{
    //Compute first and second derivative of plastic potential
    auto dGDInv = computeDGDInv(sig, rho, tempKappa);
    auto dDGDInvDKappa = computeDDGDInvDKappa(sig, rho, tempKappa);

    double equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv [ 0 ], 2.) + pow(dGDInv [ 1 ], 2.) );

    double ductilityMeasure = computeDuctilityMeasure(sig, rho, theta);
    //Compute dEquivalentDGDStressDKappa
    double dEquivalentDGDStressDKappa =
        ( 2. / 3. * dGDInv [ 0 ] * dDGDInvDKappa [ 0 ] + 2. * dGDInv [ 1 ] * dDGDInvDKappa [ 1 ] ) / ( 2. * equivalentDGDStress );

    return dEquivalentDGDStressDKappa / ductilityMeasure;
}

FloatArrayF< 2 >
ConcreteDPM2::computeDDuctilityMeasureDInv(double sig,
                                           double rho,
                                           double theta,
                                           double tempKappa) const
{
    double thetaConst = pow(2. * cos(theta), 2.);
    double x = ( -( sig + fc / 3. ) ) / fc;

    if ( x < 0. ) {
        double dXDSig = -1. / fc;
        /* Introduce exponential help function which results in a
         * smooth transition. */
        double EHard = BHard - DHard;
        double FHard = ( BHard - DHard ) * CHard / ( AHard - BHard );

        double dDuctilityMeasureDX = EHard / FHard * exp(x / FHard) / thetaConst;
        return {
            dDuctilityMeasureDX *dXDSig, 0.
        };
    } else {
        double dXDSig = -1. / fc;
        double dDuctilityMeasureDX = -( BHard - AHard ) / ( CHard ) / thetaConst * exp(-x / ( CHard ) );
        return {
            dDuctilityMeasureDX *dXDSig, 0.
        };
    }
}

FloatArrayF< 2 >
ConcreteDPM2::computeDGDInv(double sig,
                            double rho,
                            double tempKappa) const
{
    //compute yieldHard and yieldSoft
    double yieldHardOne = computeHardeningOne(tempKappa);
    double yieldHardTwo = computeHardeningTwo(tempKappa);

    //Compute dilation parameter
    double AGParam = this->ft * yieldHardTwo * 3 / this->fc + m / 2;
    double BGParam =
        yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) /
        ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

    double R = ( sig - ft / 3. * yieldHardTwo ) / fc / BGParam;
    double mQ = AGParam * exp(R);

    double Bl = sig / fc + rho / ( fc * sqrt(6.) );

    double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fc;

    double dgdsig = 4. * ( 1. - yieldHardOne ) / fc * Al * Bl + pow(yieldHardOne, 2.) * mQ / fc;

    double dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +
                    m * pow(yieldHardOne, 2.) / ( sqrt(6.) * fc );

    return {
        dgdsig, dgdrho
    };
}

FloatArrayF< 6 >
ConcreteDPM2::computeDFDStress(const FloatArrayF< 6 > &stress,
                               double tempKappa) const
{
    double sig, rho, theta;
    computeCoordinates(stress, sig, rho, theta);
    auto dFDInv = computeDFDInv(sig, rho, theta, tempKappa);

    double dRDCosTheta = computeDRDCosTheta(theta, this->ecc);

    double yieldHardOne = computeHardeningOne(tempKappa);
    double yieldHardTwo = computeHardeningTwo(tempKappa);

    //Compute dFDCosTheta. This was not needed for the stress return, but now for the tangent stiffness
    auto dFDCosTheta = dRDCosTheta * pow(yieldHardOne, 2.) * yieldHardTwo * m * rho / ( sqrt(6.) * fc );

    auto dSigDStress = computeDSigDStress();
    auto dRhoDStress = computeDRhoDStress(stress);
    auto dCosThetaDStress = computeDCosThetaDStress(stress);

    dSigDStress *= dFDInv.at(1);
    dRhoDStress *= dFDInv.at(2);
    dCosThetaDStress *= dFDCosTheta;

    auto dFDStress = dSigDStress + dRhoDStress + dCosThetaDStress;

    return dFDStress;
}

FloatArrayF< 6 >
ConcreteDPM2::computeDGDStress(const FloatArrayF< 6 > &stress, const double tempKappa) const
{
    auto tmp = computeDeviatoricVolumetricSplit(stress);
    auto deviatoricStress = tmp.first;
    double sig = tmp.second;
    double rho = computeSecondCoordinate(deviatoricStress);

    //compute dGDSig*dSigDStress + dGDRho*dRhoDStress
    auto dGDInv = computeDGDInv(sig, rho, tempKappa);
    auto dSigDStress = computeDSigDStress();
    auto dRhoDStress = computeDRhoDStress(stress);

    dSigDStress *= dGDInv.at(1);
    dRhoDStress *= dGDInv.at(2);
    dSigDStress += dRhoDStress;

    return dSigDStress;
}



//Debug
FloatMatrixF< 8, 8 >
ConcreteDPM2::computeFullJacobian(const FloatArrayF< 6 > &stress,
                                  const double deltaLambda,
                                  GaussPoint *gp,
                                  TimeStep *atTime,
                                  const double tempKappa) const
{
    FloatMatrixF< 8, 8 >jacobian;
    auto dFDStress = computeDFDStress(stress, tempKappa);
    auto dGDStress = computeDGDStress(stress, tempKappa);
    auto dDGDDStress = computeDDGDDStress(stress, tempKappa);
    auto dDGDStressDKappa = computeDDGDStressDKappa(stress, tempKappa);
    auto dDKappaDDeltaLambdaDStress = computeDDKappaDDeltaLambdaDStress(stress, tempKappa);

    double sig, rho, theta;
    computeCoordinates(stress, sig, rho, theta);

    double dFDKappa = computeDFDKappa(sig, rho, theta, tempKappa);
    double dKappaDDeltaLambda  = computeDKappaDDeltaLambda(sig, rho, theta, tempKappa);

    double dDKappaDDeltaLambdaDKappa = computeDDKappaDDeltaLambdaDKappa(sig, rho, theta, tempKappa);

    auto d = this->linearElasticMaterial.give3dMaterialStiffnessMatrix(ElasticStiffness, gp, atTime);

    auto helpA = dot(d, dDGDDStress);

    //Assign jacobian
    for ( int i = 0; i < 6; i++ ) {
        for ( int j = 0; j < 6; j++ ) {
            if ( i == j ) {
                jacobian.at(i + 1, j + 1) = 1. + deltaLambda * helpA.at(i + 1, j + 1);
            } else {
                jacobian.at(i + 1, j + 1) = deltaLambda * helpA.at(i + 1, j + 1);
            }
        }
    }

    FloatArrayF< 6 >helpB;
    helpB = dot(d, dDGDStressDKappa);

    for ( int i = 0; i < 6; i++ ) {
        jacobian.at(i + 1, 7) = deltaLambda * helpB.at(i + 1);
    }

    helpB = dot(d, dGDStress);

    for ( int i = 0; i < 6; i++ ) {
        jacobian.at(i + 1, 8) = helpB.at(i + 1);
    }

    for ( int i = 0; i < 6; i++ ) {
        jacobian.at(7, i + 1) = deltaLambda * dDKappaDDeltaLambdaDStress.at(i + 1);
    }


    jacobian.at(7, 7) = deltaLambda * dDKappaDDeltaLambdaDKappa - 1.;
    jacobian.at(7, 8) = dKappaDDeltaLambda;

    for ( int i = 0; i < 6; i++ ) {
        jacobian.at(8, i + 1) = dFDStress.at(i + 1);
    }

    jacobian.at(8, 7) = dFDKappa;
    jacobian.at(8, 8) = 0.;

    return jacobian;
}


FloatMatrixF< 6, 6 >
ConcreteDPM2::computeDDGDDStress(const FloatArrayF< 6 > &stress,
                                 const double tempKappa) const
{
    FloatMatrixF< 6, 6 >answer;

    auto tmp = computeDeviatoricVolumetricSplit(stress);
    auto deviatoricStress = tmp.first;
    double sig = tmp.second;
    double rho = computeSecondCoordinate(deviatoricStress);

    //compute deriviates with respect to invariants
    auto dDGDDInv = computeDDGDDInv(sig, rho, tempKappa);
    auto dGDInv = computeDGDInv(sig, rho, tempKappa);

    //Compute derivatives of invariants with respect to stress
    auto dRhoDStress =  computeDRhoDStress(stress);

    auto dDRhoDDStress = computeDDRhoDDStress(stress);

    auto dSigDStress =  computeDSigDStress();

    //Assemble terms

    //Compute (dDGDDSig*dSigDStress + dDGDSigDRho*dRhoDStress)*dSigDStress
    FloatArrayF< 6 >temp1;
    temp1 = dSigDStress;
    temp1 *= dDGDDInv.at(1, 1);

    FloatArrayF< 6 >temp2;
    temp2 = dRhoDStress;
    temp2 *= dDGDDInv.at(1, 2);

    temp1 += temp2;

    FloatMatrixF< 6, 6 >helpA;
    helpA = dyad(temp1, dSigDStress);

    //Compute (dDGDDRho*dRhoDStress + dDGDRhoDSig*dSigDStress)*dRhoDstress
    temp1 = dRhoDStress;
    temp1 *= dDGDDInv.at(2, 2);

    temp2 = dSigDStress;
    temp2 *= dDGDDInv.at(2, 1);

    temp1 += temp2;

    FloatMatrixF< 6, 6 >helpB;
    helpB = dyad(temp1, dRhoDStress);

    //compute dGDRho * dDRhoDDStress
    FloatMatrixF< 6, 6 >helpC = dDRhoDDStress;
    helpC *= dGDInv.at(2);

    //The term dGDSig*dDSigDDStress is zero

    //sum up all parts
    answer = helpA;
    answer += helpB;
    answer += helpC;

    return answer;
}

FloatArrayF< 2 >
ConcreteDPM2::computeDDGDInvDKappa(double sig,
                                   double rho,
                                   double tempKappa) const
{
    //Compute dilation parameter

    //compute yieldHard and yieldSoft
    double yieldHardOne = computeHardeningOne(tempKappa);
    double yieldHardTwo = computeHardeningTwo(tempKappa);

    double dYieldHardOneDKappa = computeHardeningOnePrime(tempKappa);
    double dYieldHardTwoDKappa = computeHardeningTwoPrime(tempKappa);

    //Compute dilation parameter
    double AGParam = this->ft * yieldHardTwo * 3 / this->fc + m / 2;
    double BGParam =
        yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) /
        ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

    double R = ( sig - ft / 3. * yieldHardTwo ) / ( fc * BGParam );
    double mQ = AGParam * exp(R);

    //Compute the derivative of mQ with respect to kappa

    //Derivative of AGParam
    double dAGParamDKappa = dYieldHardTwoDKappa * 3. * this->ft / this->fc;

    //Derivative of BGParam
    double BGParamTop = yieldHardTwo / 3. * ( 1. + this->ft / this->fc );
    double BGParamBottom = ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

    double dBGParamTopDKappa = dYieldHardTwoDKappa / 3. * ( 1. + this->ft / this->fc );
    double dBGParamBottomDKappa = 1. / AGParam * dAGParamDKappa - 3. * dYieldHardTwoDKappa / ( 3 * yieldHardTwo + m / 2. );
    double dBGParamDKappa = ( dBGParamTopDKappa * BGParamBottom - BGParamTop * dBGParamBottomDKappa ) / pow(BGParamBottom, 2.);

    //Derivative of R
    double RTop = ( sig - ft / 3. * yieldHardTwo );
    double RBottom = fc * BGParam;
    double dRTopDKappa = -this->ft / 3. * dYieldHardTwoDKappa;
    double dRBottomDKappa = this->fc * dBGParamDKappa;
    double dRDKappa = ( dRTopDKappa * RBottom - RTop * dRBottomDKappa ) / pow(RBottom, 2.);

    double dMQDKappa = dAGParamDKappa * exp(R) + AGParam * dRDKappa * exp(R);

    double Bl = sig / fc + rho / ( fc * sqrt(6.) );

    double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fc;

    double dAlDYieldHard = -pow(Bl, 2.);

    const double dDGDSigDKappa =
        ( -4. * Al * Bl / fc + 4. * ( 1 - yieldHardOne ) / fc * dAlDYieldHard * Bl ) * dYieldHardOneDKappa +
        dYieldHardOneDKappa * 2 * yieldHardOne * mQ / fc + pow(yieldHardOne, 2.) * dMQDKappa / fc;

    //dDGDSigDKappa = 0.;

    const double dDGDRhoDKappa =
        ( dAlDYieldHard / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) -
          4. * Al / ( sqrt(6.) * fc ) * Bl + m / ( sqrt(6.) * fc ) * 2 * yieldHardOne ) * dYieldHardOneDKappa;

    return {
        dDGDSigDKappa, dDGDRhoDKappa
    };
}

FloatMatrixF< 2, 2 >
ConcreteDPM2::computeDDGDDInv(double sig,
                              double rho,
                              double tempKappa) const
{
    //compute yieldHardOne and yieldSoft
    double yieldHardOne = computeHardeningOne(tempKappa);
    double yieldHardTwo = computeHardeningTwo(tempKappa);

    //CoQpute dilation parameter
    double AGParam = this->ft * yieldHardTwo * 3 / this->fc + m / 2;
    double BGParam =
        yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) /
        ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

    double R = ( sig - ft / 3. * yieldHardTwo ) / fc / BGParam;

    double dMQDSig = AGParam / ( BGParam * fc ) * exp(R);

    //compute help parameter Al and Bl and the corresponding derivatives
    double Bl = sig / fc + rho / ( fc * sqrt(6.) );

    double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) +
                sqrt(3. / 2.) * rho / fc;

    double dAlDSig = 2. * ( 1. - yieldHardOne ) * Bl / fc;
    double dBlDSig = 1. / fc;

    double dAlDRho = 2. * ( 1. - yieldHardOne ) * Bl / ( fc * sqrt(6.) ) + sqrt(3. / 2.) / fc;
    double dBlDRho = 1. / ( fc * sqrt(6.) );

    //compute second derivatives of g
    double ddgddSig = 4. * ( 1. - yieldHardOne ) / fc * ( dAlDSig * Bl + Al * dBlDSig ) +
                      pow(yieldHardOne, 2.) * dMQDSig / fc;
    // ddgddSig = 0.;

    double ddgddRho = dAlDRho / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +
                      Al * dBlDRho * 4. * ( 1. - yieldHardOne ) / ( sqrt(6.) * fc );

    double ddgdSigdRho = 4. * ( 1. - yieldHardOne ) / fc * ( dAlDRho * Bl + Al * dBlDRho );
    // ddgdSigdRho = 0.;

    double ddgdRhodSig = dAlDSig / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. )
                         + Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * dBlDSig );
    // ddgdRhodSig = 0.;

    FloatMatrixF< 2, 2 >answer;
    answer.at(1, 1) = ddgddSig;
    answer.at(1, 2) = ddgdSigdRho;
    answer.at(2, 1) = ddgdRhodSig;
    answer.at(2, 2) = ddgddRho;
    return answer;
}



double
ConcreteDPM2::computeAlpha(FloatArrayF< 6 > &effectiveStressTension,
                           FloatArrayF< 6 > &effectiveStressCompression,
                           const FloatArrayF< 6 > &effectiveStress) const
{
    auto tmp = StructuralMaterial::computePrincipalValDir(from_voigt_stress(effectiveStress) );
    auto principalStress = tmp.first;
    auto stressPrincipalDir = tmp.second;

    //Split the principal values in a tension and a compression part
    FloatArrayF< 6 >principalStressTension;
    FloatArrayF< 6 >principalStressCompression;

    for ( int i = 1; i <= principalStress.giveSize(); i++ ) {
        if ( principalStress.at(i) >= 0 ) {
            principalStressTension.at(i) = principalStress.at(i);
            principalStressCompression.at(i) = 0.;
        } else {
            principalStressTension.at(i) = 0.;
            principalStressCompression.at(i) = principalStress.at(i);
        }
    }

    //Transform the tension and compression principal stresses back to the original coordinate system

    //Take care of type of stress state for tension
    effectiveStressTension = transformStressVectorTo(stressPrincipalDir, principalStressTension, 1);

    //Take care of type of stress state for compression
    effectiveStressCompression = transformStressVectorTo(stressPrincipalDir, principalStressCompression, 1);

    //Determine the two factors from the stress

    double squareNormOfPrincipalStress = norm_squared(principalStress);

    double alphaTension = 0.;

    if ( squareNormOfPrincipalStress > 0 ) {
        for ( int i = 1; i <= principalStress.giveSize(); i++ ) {
            alphaTension += principalStressTension.at(i) *
                            ( principalStressTension.at(i) + principalStressCompression.at(i) ) / squareNormOfPrincipalStress;
        }
    }

    return 1. - alphaTension;
}


double
ConcreteDPM2::computeHardeningOne(double kappa) const
{
    if ( kappa <= 0. ) {
        return yieldHardInitial;
    } else if ( kappa > 0. && kappa < 1. ) {
        return
            ( 1. - yieldHardInitial - yieldHardPrimePeak ) * pow(kappa, 3.)
            - ( 3. * ( 1. - yieldHardInitial ) - 3. * yieldHardPrimePeak ) * pow(kappa, 2.)
            + ( 3. * ( 1. - yieldHardInitial ) - 2. * yieldHardPrimePeak ) * kappa
            + yieldHardInitial;
    } else {
        return 1.;
    }
}


double
ConcreteDPM2::computeHardeningOnePrime(double kappa) const
{
    if ( kappa <= 0. ) {
        return 3. * ( 1 - yieldHardInitial ) - 2. * yieldHardPrimePeak;
    } else if ( kappa >= 0. && kappa < 1. ) {
        return
            3. * ( 1. - yieldHardInitial - yieldHardPrimePeak ) * pow(kappa, 2.)
            - 2. * ( 3. * ( 1. - yieldHardInitial ) - 3. * yieldHardPrimePeak ) * kappa
            + ( 3. * ( 1. - yieldHardInitial ) - 2. * yieldHardPrimePeak );
    } else {
        return 0.;
    }
}


double
ConcreteDPM2::computeHardeningTwo(double kappa) const
{
    if ( kappa <= 0. ) {
        return 1.;
    } else if ( kappa > 0. && kappa < 1. ) {
        return 1.;
    } else {
        return 1. + ( kappa - 1. ) * yieldHardPrimePeak;
    }
}


double
ConcreteDPM2::computeHardeningTwoPrime(double kappa) const
{
    if ( kappa <= 0. ) {
        return 0.;
    } else if ( kappa >= 0. && kappa < 1. ) {
        return 0.;
    } else {
        return yieldHardPrimePeak;
    }
}

FloatMatrixF< 6, 6 >
ConcreteDPM2::give3dMaterialStiffnessMatrix(MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep) const
{
    auto status = static_cast< ConcreteDPM2Status * >( giveStatus(gp) );
    if ( mode == SecantStiffness ) {
        return this->compute3dSecantStiffness(gp, tStep);
    } else if ( mode == TangentStiffness ) {
        const int stateFlag = status->giveTempStateFlag();
        const int returnType = status->giveTempReturnType();
        if ( ( stateFlag == ConcreteDPM2Status::ConcreteDPM2_PlasticDamage || stateFlag == ConcreteDPM2Status::ConcreteDPM2_Plastic ) && returnType == ConcreteDPM2Status::RT_Regular ) {
            return this->compute3dTangentStiffness(gp, tStep);
        } else {
            return this->compute3dSecantStiffness(gp, tStep);
        }
    }
    //Elastic case
    return this->linearElasticMaterial.give3dMaterialStiffnessMatrix(mode, gp, tStep);
}


FloatMatrixF< 6, 6 >
ConcreteDPM2::compute3dSecantStiffness(GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< ConcreteDPM2Status * >( giveStatus(gp) );

    //  Damage parameters
    double omegaTension = min(status->giveTempDamageTension(), 0.999999);
    double omegaCompression = min(status->giveTempDamageCompression(), 0.999999);

    double alpha = status->giveTempAlpha();

    auto d = this->linearElasticMaterial.give3dMaterialStiffnessMatrix(ElasticStiffness, gp, tStep);

    if ( damageFlag == 2 ) {
        d *= ( 1. - omegaTension );
    } else if ( damageFlag == 3 ) {
        d *= ( 1. - ( 1. - alpha ) * omegaTension ) * ( 1. - alpha * omegaCompression );
    }

    return d;
}

FloatMatrixF< 6, 6 >
ConcreteDPM2::compute3dTangentStiffness(GaussPoint *gp, TimeStep *tStep) const
{
    FloatMatrixF< 6, 6 >answer;

    FloatArrayF< 6 >effectiveStress;
    auto status = static_cast< ConcreteDPM2Status * >( giveStatus(gp) );
    effectiveStress = status->giveTempEffectiveStress();

    double tempKappa = status->giveTempKappaP();
    double deltaLambda = status->giveDeltaLambda();

    //Computes only the plastic part of the tangent stiffness
    auto d = this->linearElasticMaterial.give3dMaterialStiffnessMatrix(ElasticStiffness, gp, tStep);

    // Debug
    FloatMatrixF< 8, 8 >fullJacobian = computeFullJacobian(effectiveStress, deltaLambda, gp, tStep, tempKappa);

    FloatMatrixF< 8, 8 >invFullJacobian = inv(fullJacobian);

    FloatMatrixF< 6, 6 >help;
    for ( int i = 1; i <= 6; i++ ) {
        for ( int j = 1; j <= 6; j++ ) {
            help.at(i, j) = invFullJacobian.at(i, j);
        }
    }

    answer = dot(help, d);

    //  Damage parameters
    double omegaTension = min(status->giveTempDamageTension(), 0.999999);
    double omegaCompression = min(status->giveTempDamageCompression(), 0.999999);
    double alpha = status->giveTempAlpha();

    if ( damageFlag == 2 ) {
        answer *= ( 1. - ( 1. - alpha ) * omegaTension ) * ( 1. - alpha * omegaCompression );
    } else if ( damageFlag == 3 ) {
        answer *= ( 1. - omegaTension );
    }

    return answer;
}

void
ConcreteDPM2::computeCoordinates(const FloatArrayF< 6 > &stress, double &sigNew, double &rhoNew, double &thetaNew) const
{
    auto tmp = computeDeviatoricVolumetricSplit(stress);
    auto deviatoricStress = tmp.first;
    sigNew = tmp.second;
    rhoNew = computeSecondCoordinate(deviatoricStress);
    thetaNew = computeThirdCoordinate(deviatoricStress);
}

void
ConcreteDPM2::assignStateFlag(GaussPoint *gp) const
{
    auto status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    //Get kappaD from status to define state later on
    double damageTension = status->giveDamageTension();
    double damageCompression = status->giveDamageCompression();
    double tempDamageTension = status->giveTempDamageTension();
    double tempDamageCompression = status->giveTempDamageCompression();
    double kappaP = status->giveKappaP();
    double tempKappaP = status->giveTempKappaP();

    if ( tempKappaP > kappaP ) {
        if ( tempDamageTension > damageTension ||  tempDamageTension == 1. ||
             tempDamageCompression > damageCompression || tempDamageCompression == 1. ) {
            status->letTempStateFlagBe(ConcreteDPM2Status::ConcreteDPM2_PlasticDamage);
        } else {
            status->letTempStateFlagBe(ConcreteDPM2Status::ConcreteDPM2_Plastic);
        }
    } else {
        const int state_flag = status->giveStateFlag();
        if ( tempDamageTension > damageTension || tempDamageTension == 1. ||
             tempDamageCompression > damageCompression || tempDamageCompression == 1. ) {
            status->letTempStateFlagBe(ConcreteDPM2Status::ConcreteDPM2_Damage);
        } else {
            if ( state_flag == ConcreteDPM2Status::ConcreteDPM2_Elastic ) {
                status->letTempStateFlagBe(ConcreteDPM2Status::ConcreteDPM2_Elastic);
            } else {
                status->letTempStateFlagBe(ConcreteDPM2Status::ConcreteDPM2_Unloading);
            }
        }
    }
}


double
ConcreteDPM2::computeDRDCosTheta(const double theta, const double ecc) const
{
    //    double thetaCrit = M_PI/3.*0.999;
    double ACostheta = 4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta) +
                       ( 2. * ecc - 1. ) * ( 2. * ecc - 1. );
    double BCostheta = 2. * ( 1. - ecc * ecc ) * cos(theta) +
                       ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta)
                                                + 5. * ecc * ecc - 4. * ecc);
    double A1Costheta = 8. * ( 1. - pow(ecc, 2.) ) * cos(theta);
    double B1Costheta = 2. * ( 1. - pow(ecc, 2.) ) +
                        4. * ( 2. * ecc - 1. ) * ( 1. - pow(ecc, 2.) ) * cos(theta) /
                        sqrt(4. * ( 1. - pow(ecc, 2.) ) * pow(cos(theta), 2.) +
                             5. * pow(ecc, 2.) - 4. * ecc);
    double dRDCostheta = A1Costheta / BCostheta - ACostheta / pow(BCostheta, 2.) * B1Costheta;
    return dRDCostheta;
}

double
ConcreteDPM2::computeDDRDDCosTheta(const double theta, const double ecc) const
{
    //  double thetaCrit = M_PI/3.*0.999;
    // if (thetaCrit<theta){
    //  printf("theta = %e\n", theta);
    //  return 0;
    // }

    //compute the derivative of the  rFunction with respect to costheta
    double atheta = 4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta) +
                    ( 2. * ecc - 1. ) * ( 2. * ecc - 1. );
    double btheta = 2. * ( 1. - ecc * ecc ) * cos(theta) +
                    ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta)
                                             + 5. * ecc * ecc - 4. * ecc);
    double a1theta = 8. * ( 1. - ecc * ecc ) * cos(theta);
    double b1theta = 2. * ( 1. - ecc * ecc ) +
                     4. * ( 2. * ecc - 1. ) * ( 1. - ecc * ecc ) * cos(theta) /
                     sqrt(4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta) +
                          5. * ecc * ecc - 4. * ecc);

    //compute the second derivative of rFunction with respect to theta
    double a2theta = 8. * ( 1. - pow(ecc, 2.) );
    double Ntheta = 4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta) +
                    5. * ecc * ecc - 4. * ecc;
    printf( "cos(theta) = %e\n", cos(theta) );
    double b2theta =  4. * ( 2. * ecc - 1. ) * ( 1. - ecc * ecc ) / sqrt(Ntheta) *
                     ( 1. - 4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta) / Ntheta );
    double ddrddcostheta = a2theta / btheta - 2. * a1theta * b1theta / ( btheta * btheta ) -
                           atheta * b2theta / ( btheta * btheta ) +
                           2. * atheta * b1theta * b1theta / ( btheta * btheta * btheta );

    return ddrddcostheta;
}


FloatArrayF< 6 >
ConcreteDPM2::computeDCosThetaDStress(const FloatArrayF< 6 > &stress) const
{
    //compute volumetric deviatoric split
    auto deviatoricStress = computeDeviator(stress);

    //compute the coordinates
    double rho = computeSecondCoordinate(deviatoricStress);

    //compute principal stresses and directions
    //auto [principalDeviatoricStress, principalDir] = computePrincipalValDir(from_voigt_stress(deviatoricStress)); // c++17
    auto tmp = computePrincipalValDir(from_voigt_stress(deviatoricStress) );
    auto principalDeviatoricStress = tmp.first;
    auto principalDir = tmp.second;

    //compute the derivative of s1 with respect to the cartesian stress
    FloatArrayF< 6 >ds1DStress;
    ds1DStress.at(1) = principalDir.at(1, 1) * principalDir.at(1, 1) - 1. / 3.;
    ds1DStress.at(2) = principalDir.at(2, 1) * principalDir.at(2, 1) - 1. / 3.;
    ds1DStress.at(3) = principalDir.at(3, 1) * principalDir.at(3, 1) - 1. / 3.;
    ds1DStress.at(4) = 2. * principalDir.at(2, 1) * principalDir.at(3, 1);
    ds1DStress.at(5) = 2. * principalDir.at(3, 1) * principalDir(1, 1);
    ds1DStress.at(6) = 2. * principalDir.at(1, 1) * principalDir.at(2, 1);

    //compute dCosThetaDStress
    auto dCosThetaDStress = ds1DStress * ( sqrt(3. / 2.) * rho / pow(rho, 2.) ) +
                            computeDRhoDStress(stress) * ( -sqrt(3. / 2.) * principalDeviatoricStress [ 0 ] / pow(rho, 2.) );
    return dCosThetaDStress;
}


FloatMatrixF< 6, 6 >
ConcreteDPM2::computeDDCosThetaDDStress(const FloatArrayF< 6 > &stress) const
{
    //compute volumetric deviatoric split
    auto deviatoricStress = computeDeviator(stress);

    //compute the coordinate
    double rho = computeSecondCoordinate(deviatoricStress);

    //compute principal stresses and directions
    //auto [principalDeviatoricStress, principalDir] = computePrincipalValDir(from_voigt_stress(deviatoricStress)); // c++17
    auto tmp = computePrincipalValDir(from_voigt_stress(deviatoricStress) );
    auto principalDeviatoricStress = tmp.first;
    auto principalDir = tmp.second;


    //compute the derivative of s1 with respect to the cartesian stress
    FloatArrayF< 6 >dS1DStress;
    dS1DStress.at(1) = principalDir.at(1, 1) * principalDir.at(1, 1) - 1. / 3.;
    dS1DStress.at(2) = principalDir.at(2, 1) * principalDir.at(2, 1) - 1. / 3.;
    dS1DStress.at(3) = principalDir.at(3, 1) * principalDir.at(3, 1) - 1. / 3.;
    dS1DStress.at(4) = 2. * principalDir.at(2, 1) * principalDir.at(3, 1);
    dS1DStress.at(5) = 2. * principalDir.at(3, 1) * principalDir.at(1, 1);
    dS1DStress.at(6) = 2. * principalDir.at(1, 1) * principalDir.at(2, 1);

    //compute the second derivative of costheta
    FloatMatrixF< 6, 6 >dDCosThetaDDStress;
    //compute the product of dS1DStress and DRhoDStress
    auto dRhoDStress = computeDRhoDStress(stress);
    FloatMatrixF< 6, 6 >dS1DRho;
    for ( int v = 0; v < 6; v++ ) {
        for ( int w = 0; w < 6; w++ ) {
            dS1DRho.at(v + 1, w + 1) = dS1DStress.at(v + 1) * dRhoDStress.at(w + 1);
        }
    }

    //compute the square of dRhoDStress and DRhoDStress
    FloatMatrixF< 6, 6 >dRhoDRho;
    for ( int v = 0; v < 6; v++ ) {
        for ( int w = 0; w < 6; w++ ) {
            dRhoDRho.at(v + 1, w + 1) = dRhoDStress.at(v + 1) * dRhoDStress.at(w + 1);
        }
    }

    //compute the product of dRhoDStress and DS1DStress
    FloatMatrixF< 6, 6 >dRhoDS1;
    for ( int v = 0; v < 6; v++ ) {
        for ( int w = 0; w < 6; w++ ) {
            dRhoDS1.at(v + 1, w + 1) = dRhoDStress.at(v + 1) * dS1DStress.at(w + 1);
        }
    }

    auto dDRhoDDStress = computeDDRhoDDStress(stress);
    dDCosThetaDDStress = dDRhoDDStress;
    dDCosThetaDDStress *= ( -sqrt(3. / 2.) * principalDeviatoricStress.at(1) / pow(rho, 2.) );

    auto helpB1 = dS1DRho;
    helpB1 *= ( -sqrt(3. / 2.) / pow(rho, 2.) );

    auto helpC1 = dRhoDRho;
    helpC1 *= ( sqrt(6.) * principalDeviatoricStress.at(1) / pow(rho, 3.) );

    auto helpD1 = dRhoDS1;
    helpD1 *= ( -sqrt(3. / 2.) / pow(rho, 2.) );
    dDCosThetaDDStress += helpB1;
    dDCosThetaDDStress += helpC1;
    dDCosThetaDDStress += helpD1;

    return dDCosThetaDDStress;
}


FloatArrayF< 6 >
ConcreteDPM2::computeDRhoDStress(const FloatArrayF< 6 > &stress) const
{
    //compute volumetric deviatoric split
    auto deviatoricStress = computeDeviator(stress);
    double rho = computeSecondCoordinate(deviatoricStress);

    //compute the derivative of J2 with respect to the stress
    auto dJ2DStress = deviatoricStress;
    for ( int i = 3; i < 6; i++ ) {
        dJ2DStress.at(i + 1) = deviatoricStress.at(i + 1) * 2.0;
    }

    //compute the derivative of rho with respect to stress
    return dJ2DStress * ( 1. / rho );
}

FloatArrayF< 6 >
ConcreteDPM2::computeDSigDStress() const
{
    return {
        1. / 3., 1. / 3., 1. / 3., 0., 0., 0.
    };
}


FloatMatrixF< 6, 6 >
ConcreteDPM2::computeDDRhoDDStress(const FloatArrayF< 6 > &stress) const
{
    //compute volumetric deviatoric split
    auto deviatoricStress = computeDeviator(stress);
    double rho = computeSecondCoordinate(deviatoricStress);

    //compute first dericative of J2
    auto dJ2dstress = deviatoricStress;
    for ( int i = 3; i < 6; i++ ) {
        dJ2dstress.at(i + 1) = deviatoricStress.at(i + 1) * 2.;
    }

    //compute second derivative of J2
    FloatMatrixF< 6, 6 >ddJ2ddstress;
    for ( int i = 0; i < 6; i++ ) {
        if ( i < 3 ) {
            ddJ2ddstress.at(i + 1, i + 1) = 2. / 3.;
        }

        if ( i > 2 ) {
            ddJ2ddstress.at(i + 1, i + 1) = 2.;
        }
    }

    ddJ2ddstress.at(1, 2) = -1. / 3.;
    ddJ2ddstress.at(1, 3) = -1. / 3.;
    ddJ2ddstress.at(2, 1) = -1. / 3.;
    ddJ2ddstress.at(2, 3) = -1. / 3.;
    ddJ2ddstress.at(3, 1) = -1. / 3.;
    ddJ2ddstress.at(3, 2) = -1. / 3.;

    return ddJ2ddstress * ( 1. / rho ) + dyad(dJ2dstress, dJ2dstress) * ( -1. / ( rho * rho * rho ) );
}

int
ConcreteDPM2::giveIPValue(FloatArray &answer,
                          GaussPoint *gp,
                          InternalStateType type,
                          TimeStep *tStep)
{
    auto status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    switch ( type ) {
    case IST_PlasticStrainTensor:
        answer = status->givePlasticStrain();
        return 1;

    case IST_DamageTensor:
        answer.resize(6);
        answer.zero();
        answer.at(1) = status->giveDamageTension();
        answer.at(2) = status->giveDamageCompression();
        answer.at(3) = status->giveDissWork();
        return 1;

#ifdef keep_track_of_dissipated_energy
    case IST_StressWorkDensity:
        answer.resize(1);
        answer.at(1) = status->giveStressWork();
        return 1;

    case IST_DissWorkDensity:
        answer.resize(1);
        answer.at(1) = status->giveDissWork();
        return 1;

    case IST_FreeEnergyDensity:
        answer.resize(1);
        answer.at(1) = status->giveStressWork() - status->giveDissWork();
        return 1;

#endif
    default:
        return StructuralMaterial::giveIPValue(answer, gp, type, tStep);
    }
}

MaterialStatus *
ConcreteDPM2::CreateStatus(GaussPoint *gp) const
{
    return new  ConcreteDPM2Status(gp);
}
} //end of namespace
