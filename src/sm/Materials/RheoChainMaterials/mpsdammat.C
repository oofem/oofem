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

#include "mpsdammat.h"
#include "mathfem.h"
#include "gausspoint.h"
#include "timestep.h"
#include "stressvector.h"
#include "engngm.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(MPSDamMaterial);

/***************************************************************/
/**************     MPSDamMaterialStatus     *******************/
/***************************************************************/



MPSDamMaterialStatus :: MPSDamMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits) :
    MPSMaterialStatus(n, d, g, nunits), effectiveStressVector(), tempEffectiveStressVector()
{
    kappa = tempKappa = 0.0;
    damage = tempDamage = 0.0;
    charLength = 0.0;
    crackVector.resize(3);
    crackVector.zero();

    var_e0 = var_gf = 0.;

    int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
    effectiveStressVector.resize(rsize);
    effectiveStressVector.zero();
    tempEffectiveStressVector = effectiveStressVector;

#ifdef supplementary_info
    crackWidth = 0.;
    residTensileStrength = 0.;
#endif
}


void
MPSDamMaterialStatus :: initTempStatus()
{
    MPSMaterialStatus :: initTempStatus();

    this->tempKappa = this->kappa;
    this->tempDamage = this->damage;

    this->tempEffectiveStressVector = this->effectiveStressVector;

    if ( !damage ) {
        var_e0 = var_gf = 0.;
    }
}

void
MPSDamMaterialStatus :: updateYourself(TimeStep *tStep)
{
    MPSMaterialStatus :: updateYourself(tStep);

    this->kappa = this->tempKappa;
    this->damage = this->tempDamage;

    this->effectiveStressVector = tempEffectiveStressVector;

    if ( !damage ) {
        var_e0 = var_gf = 0.;
    }
}

void
MPSDamMaterialStatus :: giveCrackVector(FloatArray &answer)
{
    answer = crackVector;
    answer.times(damage);
}


void
MPSDamMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    MPSMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "damage status { ");
    if ( this->kappa > 0 && this->damage <= 0 ) {
        fprintf(file, "kappa %f", this->kappa);
    } else if ( this->damage > 0.0 ) {
        fprintf( file, "kappa %f, damage %f crackVector %f %f %f", this->kappa, this->damage, this->crackVector.at(1), this->crackVector.at(2), this->crackVector.at(3) );
    }
    fprintf(file, "}\n");
}



contextIOResultType
MPSDamMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
//
{
    contextIOResultType iores;

    if ( ( iores = MPSMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(charLength) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = crackVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(var_e0) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(var_gf) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = effectiveStressVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

contextIOResultType
MPSDamMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restore the state variables from a stream
//
{
    contextIOResultType iores;
    if ( ( iores = MPSMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(kappa) ) {
        return CIO_IOERR;
    }

    if ( !stream.read(damage) ) {
        return CIO_IOERR;
    }

    if ( !stream.read(charLength) ) {
        return CIO_IOERR;
    }

    if ( ( iores = crackVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(var_e0) ) {
        return CIO_IOERR;
    }

    if ( !stream.read(var_gf) ) {
        return CIO_IOERR;
    }

    if ( ( iores = effectiveStressVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    return CIO_OK;
}

//   ***********************************************************************************
//   ***   CLASS Damage extension of MPS model for creep and shrinakge of concrete   ***
//   ***********************************************************************************


MPSDamMaterial :: MPSDamMaterial(int n, Domain *d) : MPSMaterial(n, d)

{
    maxOmega = 0.999999;

    softType = ST_Exponential_Cohesive_Crack;
    ecsMethod = ECSM_Projection;
    const_e0 = 0.;
    const_gf = 0.;
    checkSnapBack = 1; //snapback check by default
}


int
MPSDamMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    return mode == _3dMat || mode == _PlaneStress || mode == _PlaneStrain || mode == _1dMat;
}



IRResultType
MPSDamMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    MPSMaterial :: initializeFrom(ir);

    this->isotropic = false;
    if ( ir->hasField(_IFT_MPSDamMaterial_isotropic) ) {
        this->isotropic = true;
    }

    // initialize dummy elastic modulus E
    //E = 20.e9 / MPSMaterial::stiffnessFactor;
    E = 1. / MPSMaterial :: computeCreepFunction(28.01, 28);

    int damageLaw = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, damageLaw, _IFT_MPSDamMaterial_damageLaw);

    this->timeDepFracturing = false;

    if ( ir->hasField(_IFT_MPSDamMaterial_timedepfracturing) ) {
        this->timeDepFracturing = true;
        //
        IR_GIVE_FIELD(ir, fib_s, _IFT_MPSDamMaterial_fib_s);
        // the same compressive strength as for the prediction using the B3 formulas

        this->gf28 = 0.;
        this->ft28 = 0.;
        
        if (  ( ir->hasField(_IFT_MPSDamMaterial_ft28) ) && ( ir->hasField(_IFT_MPSDamMaterial_gf28) ) )  {
            
            IR_GIVE_FIELD(ir, gf28, _IFT_MPSDamMaterial_gf28);
            if (gf28 < 0.) {
              OOFEM_ERROR("Fracture energy at 28 days must be positive");
            }
            IR_GIVE_FIELD(ir, ft28, _IFT_MPSDamMaterial_ft28);
            
            if (ft28 < 0.) {
              OOFEM_ERROR("Tensile strength at 28 days must be positive");
            }
            
          } else {
            IR_GIVE_OPTIONAL_FIELD(ir, gf28, _IFT_MPSDamMaterial_gf28);
            IR_GIVE_OPTIONAL_FIELD(ir, ft28, _IFT_MPSDamMaterial_ft28);
            IR_GIVE_FIELD(ir, fib_fcm28, _IFT_MPSMaterial_fc);
          }

    } else {
        double ft;

        //applies only in this class
        switch ( damageLaw ) {
        case 0:   // exponential softening - default
            IR_GIVE_FIELD(ir, ft, _IFT_MPSDamMaterial_ft);
            this->softType = ST_Exponential_Cohesive_Crack;
            IR_GIVE_FIELD(ir, const_gf, _IFT_MPSDamMaterial_gf);
            this->const_e0 = ft / E;
            break;

        case 1:   // linear softening law
            IR_GIVE_FIELD(ir, ft, _IFT_MPSDamMaterial_ft);
            this->softType = ST_Linear_Cohesive_Crack;
            IR_GIVE_FIELD(ir, const_gf, _IFT_MPSDamMaterial_gf);
            this->const_e0 = ft / E;
            break;

        case 6:
            this->softType = ST_Disable_Damage;
            break;

        default:
            OOFEM_ERROR("Softening type number %d is unknown", damageLaw);
        }
    }

    IR_GIVE_OPTIONAL_FIELD(ir, checkSnapBack, _IFT_MPSDamMaterial_checkSnapBack);

    return IRRT_OK;
}


void
MPSDamMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep)
{
    MPSDamMaterialStatus *status = static_cast< MPSDamMaterialStatus * >( this->giveStatus(gp) );

    MaterialMode mode = gp->giveMaterialMode();

    FloatArray principalStress;
    FloatMatrix principalDir;

    StressVector tempEffectiveStress(mode);
    StressVector tempNominalStress(mode);

    double f, sigma1, kappa, tempKappa = 0.0, omega = 0.0;

    // effective stress computed by the viscoelastic MPS material
    MPSMaterial :: giveRealStressVector(tempEffectiveStress, gp, totalStrain, tStep);

    if ( !this->isActivated(tStep) ) {
        FloatArray help;
        help.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        help.zero();

        status->letTempStrainVectorBe(help);
        status->letTempStressVectorBe(help);
        status->letTempViscoelasticStressVectorBe(help);

#ifdef supplementary_info
        status->setResidualTensileStrength(0.);
        status->setCrackWidth(0.);
        status->setTempKappa(0.);
        status->setTempDamage(0.);
#endif

        return;
    }

    tempEffectiveStress.computePrincipalValDir(principalStress, principalDir);

    sigma1 = 0.;
    if ( principalStress.at(1) > 0.0 ) {
        sigma1 = principalStress.at(1);
    }

    kappa = sigma1 / this->E;

    f = kappa - status->giveKappa();

    if ( f <= 0.0 ) { // damage does not grow
        tempKappa = status->giveKappa();
        omega     = status->giveDamage();

#ifdef supplementary_info
        double residualStrength = 0.;
        double e0;

        if ( ( this->timeDepFracturing ) && ( this->givee0(gp) == 0. ) ) {
            this->initDamagedFib(gp, tStep);
        }

        e0 = this->givee0(gp);

        if ( omega == 0. ) {
            residualStrength = E * e0; // undamaged material
        } else {
            double gf, ef, wf = 0., Le;
            gf = this->givegf(gp);
            Le = status->giveCharLength();

            if ( softType == ST_Exponential_Cohesive_Crack ) { // exponential softening
                wf = gf / this->E / e0; // wf is the crack opening
            } else if ( softType == ST_Linear_Cohesive_Crack ) { // linear softening law
                wf = 2. * gf / this->E / e0; // wf is the crack opening
            } else {
                OOFEM_ERROR("Gf unsupported for softening type softType = %d", softType);
            }

            ef = wf / Le; //ef is the fracturing strain

            if ( this->softType == ST_Linear_Cohesive_Crack ) {
                residualStrength = E * e0 * ( ef - tempKappa ) / ( ef - e0 );
            } else if (  this->softType == ST_Exponential_Cohesive_Crack ) {
                residualStrength = E * e0 * exp(-1. * ( tempKappa - e0 ) / ef);
            } else {
                OOFEM_ERROR("Unknown softening type for cohesive crack model.");
            }
        }

        if ( status ) {
            status->setResidualTensileStrength(residualStrength);
        }

#endif
    } else {
        // damage grows
        tempKappa = kappa;

        FloatArray crackPlaneNormal(3);
        for ( int i = 1; i <= principalStress.giveSize(); i++ ) {
            crackPlaneNormal.at(i) = principalDir.at(i, 1);
        }
        this->initDamaged(tempKappa, crackPlaneNormal, gp, tStep);
        this->computeDamage(omega, tempKappa, gp);
    }

    answer.zero();

    if ( omega > 0. ) {
        tempNominalStress = tempEffectiveStress;

        if ( this->isotropic ) {
            // convert effective stress to nominal stress
            tempNominalStress.times(1. - omega);
            answer.add(tempNominalStress);
        } else {
            // stress transformation matrix
            FloatMatrix Tstress;

            // compute principal nominal stresses by multiplying effective stresses by damage
            for ( int i = 1; i <= principalStress.giveSize(); i++ ) {
                if ( principalStress.at(i) > 0. ) {
                    // convert principal effective stress to nominal stress
                    principalStress.at(i) *= ( 1. - omega );
                }
            }

            if ( mode == _PlaneStress ) {
                principalStress.resizeWithValues(3);
                givePlaneStressVectorTranformationMtrx(Tstress, principalDir, true);
            } else {
                principalStress.resizeWithValues(6);
                giveStressVectorTranformationMtrx(Tstress, principalDir, true);
            }

            principalStress.rotatedWith(Tstress, 'n');


            if ( mode == _PlaneStress ) { // plane stress
                answer.add(principalStress);
            } else if ( this->giveSizeOfVoigtSymVector(mode) != this->giveSizeOfVoigtSymVector(_3dMat) ) { // mode = _PlaneStrain or axial symmetry
                StressVector redFormStress(mode);
                redFormStress.convertFromFullForm(principalStress, mode);
                answer.add(redFormStress);
            } else { // 3D
                answer.add(principalStress);
            }
        }
    } else {
        answer.add(tempEffectiveStress);
    }

#ifdef supplementary_info

    if ( ( omega == 0. ) || ( sigma1 <= 0 ) ) {
        status->setCrackWidth(0.);
    } else {
        FloatArray principalStrains;
        this->computePrincipalValues(principalStrains, totalStrain, principal_strain);

        double crackWidth;
        //case when the strain localizes into narrow band and after a while the stresses relax

        double strainWithoutTemperShrink = principalStrains.at(1);
        strainWithoutTemperShrink -=  status->giveTempThermalStrain();
        strainWithoutTemperShrink -=  status->giveTempDryingShrinkageStrain();
        strainWithoutTemperShrink -=  status->giveTempAutogenousShrinkageStrain();


        crackWidth = status->giveCharLength() * omega * strainWithoutTemperShrink;

        status->setCrackWidth(crackWidth);
    }

#endif

    // update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->letTempViscoelasticStressVectorBe(tempEffectiveStress);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);
}


void
MPSDamMaterial :: initDamagedFib(GaussPoint *gp, TimeStep *tStep)
{
    MPSDamMaterialStatus *status = ( MPSDamMaterialStatus * ) this->giveStatus(gp);

    if ( status->giveDamage() == 0. ) {
        double tequiv;

        tequiv = this->computeEquivalentTime(gp, tStep, 0);

        double e0 = this->computeTensileStrength(tequiv) / this->E;
        double gf = this->computeFractureEnergy(tequiv);

        status->sete0(e0);
        status->setgf(gf);
    }
}

double
MPSDamMaterial :: givee0(GaussPoint *gp)
{
    if ( this->timeDepFracturing ) {
        MPSDamMaterialStatus *status = ( MPSDamMaterialStatus * ) this->giveStatus(gp);
        return status->givee0();
    } else {
        return this->const_e0;
    }
}

double
MPSDamMaterial :: givegf(GaussPoint *gp)
{
    if ( this->timeDepFracturing ) {
        MPSDamMaterialStatus *status = ( MPSDamMaterialStatus * ) this->giveStatus(gp);
        return status->givegf();
    } else {
        return this->const_gf;
    }
}


double
MPSDamMaterial :: computeFractureEnergy(double equivalentTime)
{

    double fractureEnergy, fractureEnergy28;
    double ftm, ftm28;

    // the fracture energy has the same time evolution as the tensile strength, 
    // the direct relation to the mean value of compressive strength according to Model Code 
    // highly overestimates the initial (early age) value of the fracture energy 
    //( fractureEnergy = 73. * fcm(t)^0.18 )

    // evolution of the mean compressive strength with respect to equivalent time/age/maturity
    // returns fcm in MPa

    /*    if (this->gf28 > 0.) {
      double fcm28mod;
      fcm28mod =  pow ( this->gf28 * MPSMaterial :: stiffnessFactor / 73. , 1. / 0.18 );
      fcm = exp( fib_s * ( 1. - sqrt(28. * MPSMaterial :: lambda0 / equivalentTime) ) ) * fcm28mod;

    } else {
      fcm = exp( fib_s * ( 1. - sqrt(28. * MPSMaterial :: lambda0 / equivalentTime) ) ) * fib_fcm28;
    }

    fractureEnergy = 73. * pow(fcm, 0.18) / MPSMaterial :: stiffnessFactor;
    */

    // 1) read or estimate the 28-day value of fracture energy

    if (this->gf28 > 0.) {
      fractureEnergy28 = this->gf28;
    } else {
      fractureEnergy28 = 73. * pow(fib_fcm28, 0.18) / MPSMaterial :: stiffnessFactor;
    }

    // 2) compute the tensile strengh according to provided equivalent time

    ftm = this->computeTensileStrength(equivalentTime);
    ftm28 = this->computeTensileStrength(28. * MPSMaterial :: lambda0);

    // 3) calculate the resulting fracture energy as gf28 * ft/ft28

    fractureEnergy = fractureEnergy28 * ftm / ftm28;

    return fractureEnergy;
}

double
MPSDamMaterial :: computeTensileStrength(double equivalentTime)
{
    double fcm, ftm;


    if (this->ft28 > 0.) {
      double fcm28mod;
      fcm28mod = pow ( this->ft28 * MPSMaterial :: stiffnessFactor / 0.3e6, 3./2. ) + 8.;
      fcm = exp( fib_s * ( 1. - sqrt(28. * MPSMaterial :: lambda0 / equivalentTime) ) ) * fcm28mod;

    } else {
      // returns fcm in MPa - formula 5.1-51, Table 5.1-9
      fcm = exp( fib_s * ( 1. - sqrt(28. * MPSMaterial :: lambda0 / equivalentTime) ) ) * fib_fcm28;
    }


    // ftm adjusted according to the stiffnessFactor (MPa by default)
    if ( fcm >= 58. ) {
      ftm = 2.12 * log ( 1. + 0.1 * fcm ) * 1.e6 / MPSMaterial :: stiffnessFactor;
    } else if ( fcm <= 20. ) {
      ftm = 0.07862 * fcm * 1.e6 / MPSMaterial :: stiffnessFactor; // 12^(2/3) * 0.3 / 20 = 0.07862
    } else {
      ftm = 0.3 * pow(fcm - 8., 2. / 3.) * 1.e6 / MPSMaterial :: stiffnessFactor; //5.1-3a
    }

    /*
    if ( fcm >= 20. ) {
        ftm = 0.3 * pow(fcm - 8., 2. / 3.) * 1.e6 / MPSMaterial :: stiffnessFactor; //5.1-3a
    } else if ( fcm < 8. ) {
        // upper formula does not hold for concretes with fcm < 8 MPa
        ftm = 0.3 * pow(fcm, 2. / 3.) * 1.e6 / MPSMaterial :: stiffnessFactor;
    } else {
        // smooth transition
        ftm = 0.3 * pow(fcm - ( 8. * ( fcm - 8. ) / ( 20. - 8. ) ), 2. / 3.) * 1.e6 / MPSMaterial :: stiffnessFactor;
        }*/

    return ftm;
}


void
MPSDamMaterial :: computeDamage(double &omega, double kappa, GaussPoint *gp)
{
    if ( this->softType == ST_Disable_Damage ) { //dummy material with no damage
        omega = 0.;
    } else {
        computeDamageForCohesiveCrack(omega, kappa, gp);
    }
}

void
MPSDamMaterial :: computeDamageForCohesiveCrack(double &omega, double kappa, GaussPoint *gp)
{
    MPSDamMaterialStatus *status = NULL;

    omega = 0.0;
    double e0 = this->givee0(gp);
    double ef = 0.;

    if ( kappa > e0 ) {
        double gf = this->givegf(gp);

        double Le;
        double wf = 0.;
        if ( softType == ST_Exponential_Cohesive_Crack ) { // exponential softening
            wf = gf / this->E / e0; // wf is the crack opening
        } else if ( softType == ST_Linear_Cohesive_Crack ) { // linear softening law
            wf = 2. * gf / this->E / e0; // wf is the crack opening
        } else {
            OOFEM_ERROR("Gf unsupported for softening type softType = %d", softType);
        }

        //    MPSDamMaterialStatus *status = static_cast< MPSDamMaterialStatus * >( this->giveStatus(gp) );
        status = static_cast< MPSDamMaterialStatus * >( this->giveStatus(gp) );

        Le = status->giveCharLength();
        ef = wf / Le; //ef is the fracturing strain
        if ( ef < e0 ) { //check that no snapback occurs
            double minGf = 0.;
            OOFEM_WARNING("ef %e < e0 %e, this leads to material snapback in element %d, characteristic length %f", ef, e0, gp->giveElement()->giveNumber(), Le);

            if ( softType == ST_Exponential_Cohesive_Crack ) { //exponential softening
                minGf = this->E * e0 * e0 * Le;
            } else if ( softType == ST_Linear_Cohesive_Crack ) { //linear softening law
                minGf = this->E * e0 * e0 * Le / 2.;
            } else {
                OOFEM_WARNING("Gf unsupported for softening type softType = %d", softType);
            }

            OOFEM_WARNING("Material number %d, decrease e0, or increase Gf from %f to Gf=%f", this->giveNumber(), gf, minGf);
            if ( checkSnapBack ) {
                OOFEM_ERROR("");
            }
        }

        if ( this->softType == ST_Linear_Cohesive_Crack ) {
            if ( kappa < ef ) {
                omega = ( ef / kappa ) * ( kappa - e0 ) / ( ef - e0 );
            } else {
                omega = 1.0; //maximum omega (maxOmega) is adjusted just for stiffness matrix in isodamagemodel.C
            }
        } else if (  this->softType == ST_Exponential_Cohesive_Crack ) {
            // exponential cohesive crack - iteration needed
            double R, Lhs, help;
            int nite = 0;
            // iteration to achieve objectivity
            // we are looking for a state in which the elastic stress is equal to
            // the stress from crack-opening relation
            // ef has now the meaning of strain
            do {
                nite++;
                help = omega * kappa / ef;
                R = ( 1. - omega ) * kappa - e0 *exp(-help); //residuum
                Lhs = kappa - e0 *exp(-help) * kappa / ef; //- dR / (d omega)
                omega += R / Lhs;
                if ( nite > 40 ) {
                    OOFEM_ERROR("algorithm not converging");
                }
            } while ( fabs(R) >= e0 * MPSDAMMAT_ITERATION_LIMIT );
        } else {
            OOFEM_ERROR("Unknown softening type for cohesive crack model.");
        }

        if ( omega > 1.0 ) {
            OOFEM_ERROR("damage parameter is %f, which is greater than 1, snap-back problems", omega);
        }

        if ( omega < 0.0 ) {
            OOFEM_WARNING("damage parameter is %f, which is smaller than 0, snap-back problems", omega);
            omega = 1.;
            if ( checkSnapBack ) {
                OOFEM_ERROR("");
            }

        }
    }


#ifdef supplementary_info
    double residualStrength = 0.;

    if ( omega == 0. ) {
        residualStrength = E * e0; // undamaged material
        status = static_cast< MPSDamMaterialStatus * >( this->giveStatus(gp) );
    } else {
        if ( this->softType == ST_Linear_Cohesive_Crack ) {
            residualStrength = E * e0 * ( ef - kappa ) / ( ef - e0 );
        } else if (  this->softType == ST_Exponential_Cohesive_Crack ) {
            residualStrength = E * e0 * exp(-1. * ( kappa - e0 ) / ef);
        } else {
            OOFEM_ERROR("Unknown softening type for cohesive crack model.");
        }
    }

    if ( status ) {
        status->setResidualTensileStrength(residualStrength);
    }

#endif
}

void
MPSDamMaterial :: initDamaged(double kappa, FloatArray &principalDirection, GaussPoint *gp, TimeStep *tStep)
{
    double le = 0.;

    MPSDamMaterialStatus *status = static_cast< MPSDamMaterialStatus * >( this->giveStatus(gp) );

    if ( this->timeDepFracturing ) {
        this->initDamagedFib(gp, tStep);
    }

    double e0 = this->givee0(gp);
    double gf = this->givegf(gp);
    double wf = 0.;

    if ( softType == ST_Disable_Damage ) {
        return;
    } else if  ( softType == ST_Exponential_Cohesive_Crack ) { // exponential softening
        wf = gf / E / e0; // wf is the crack opening
    } else if ( softType == ST_Linear_Cohesive_Crack ) { //  linear softening law
        wf = 2. * gf / E / e0; // wf is the crack opening
    } else {
        OOFEM_ERROR("Gf unsupported for softening type softType = %d", softType);
    }

    if ( ( kappa > e0 ) && ( status->giveDamage() == 0. ) ) {
        status->setCrackVector(principalDirection);

        le = gp->giveElement()->giveCharacteristicSize(gp, principalDirection, ecsMethod);
        status->setCharLength(le);

        if ( gf != 0. && e0 >= ( wf / le ) ) { // case for a given fracture energy
            OOFEM_WARNING("Fracturing strain %e is lower than the elastic strain e0=%e, possible snap-back. Element number %d, wf %e, le %e. Increase fracturing strain or decrease element size by at least %f", wf / le, e0, gp->giveElement()->giveLabel(), wf, le, e0/(wf/le) );
            if ( checkSnapBack ) {
                OOFEM_ERROR("");
            }
        }
    }
}



MaterialStatus *
MPSDamMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates a new material status corresponding to this class
 */
{
    return new MPSDamMaterialStatus(1, this->giveDomain(), gp, nUnits);
}


void
MPSDamMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                MatResponseMode mode,
                                                GaussPoint *gp,
                                                TimeStep *tStep)
{
    RheoChainMaterial :: give3dMaterialStiffnessMatrix(answer, ElasticStiffness, gp, tStep);

    if ( mode == ElasticStiffness || ( mode == SecantStiffness && !this->isotropic ) ) {
        return;
    }

    double tempDamage;
    MPSDamMaterialStatus *status = ( MPSDamMaterialStatus * ) this->giveStatus(gp);
    tempDamage = min(status->giveTempDamage(), this->maxOmega);
    answer.times(1.0 - tempDamage);
}


void
MPSDamMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep)
{
    RheoChainMaterial :: givePlaneStressStiffMtrx(answer, ElasticStiffness, gp, tStep);

    if ( mode == ElasticStiffness || ( mode == SecantStiffness && !this->isotropic ) ) {
        return;
    }

    double tempDamage;
    MPSDamMaterialStatus *status = ( MPSDamMaterialStatus * ) this->giveStatus(gp);
    tempDamage = min(status->giveTempDamage(), this->maxOmega);
    answer.times(1.0 - tempDamage);
}

void
MPSDamMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep)
{
    RheoChainMaterial :: givePlaneStrainStiffMtrx(answer, ElasticStiffness, gp, tStep);

    if ( mode == ElasticStiffness || ( mode == SecantStiffness && !this->isotropic ) ) {
        return;
    }

    double tempDamage;
    MPSDamMaterialStatus *status = ( MPSDamMaterialStatus * ) this->giveStatus(gp);
    tempDamage = min(status->giveTempDamage(), this->maxOmega);
    answer.times(1.0 - tempDamage);
}


void
MPSDamMaterial :: give1dStressStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode mode,
                                        GaussPoint *gp,
                                        TimeStep *tStep)
{
    RheoChainMaterial :: give1dStressStiffMtrx(answer, ElasticStiffness, gp, tStep);

    if ( mode == ElasticStiffness || ( mode == SecantStiffness && !this->isotropic ) ) {
        return;
    }

    double tempDamage;
    MPSDamMaterialStatus *status = ( MPSDamMaterialStatus * ) this->giveStatus(gp);
    tempDamage = min(status->giveTempDamage(), this->maxOmega);
    answer.times(1.0 - tempDamage);
}


int
MPSDamMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    MPSDamMaterialStatus *status = static_cast< MPSDamMaterialStatus * >( this->giveStatus(gp) );
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
    } else if ( type == IST_CharacteristicLength ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveCharLength();
        return 1;
    } else if ( type == IST_CrackVector ) {
        answer.resize(3);
        answer.zero();
        status->giveCrackVector(answer);
        return 1;
    } else if ( type == IST_CrackWidth ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveCrackWidth();
        return 1;
    } else if ( type == IST_ResidualTensileStrength ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) =  status->giveResidualTensileStrength();
        return 1;
    } else if ( type == IST_TensileStrength ) {
        MPSDamMaterialStatus *status = static_cast< MPSDamMaterialStatus * >( this->giveStatus(gp) );
        double tequiv = status->giveEquivalentTime();
        answer.resize(1);
        answer.zero();
            if (tequiv >= this->castingTime) {
              answer.at(1) =  this->computeTensileStrength(tequiv);
            }
        return 1;
    } else if ( type == IST_CrackIndex ) {
        //ratio of real principal stress / strength. 1 if damage already occured.
        FloatArray principalStress;
        MPSDamMaterialStatus *status = static_cast< MPSDamMaterialStatus * >( this->giveStatus(gp) );
        answer.resize(1);
        answer.zero();
        if ( status->giveDamage()>0. ){
            answer.at(1)=1.;
            return 1;
        }
        //FloatArray effectiveStress = status->giveTempViscoelasticStressVector();
        //StructuralMaterial :: computePrincipalValues(principalStress, effectiveStress, principal_stress);
        StructuralMaterial :: giveIPValue(principalStress, gp, IST_PrincipalStressTensor, tStep);
        double tequiv = status->giveEquivalentTime();
            if (tequiv >= this->castingTime) {
                double ft = this->computeTensileStrength(tequiv);
                if (ft > 1.e-20 && principalStress.at(1)>1.e-20){
                answer.at(1) = principalStress.at(1)/ft;
                } 
        }
        return 1;
    } else {
        return MPSMaterial :: giveIPValue(answer, gp, type, tStep);
    }

    return 1; // to make the compiler happy
}
}  // end namespace oofem
