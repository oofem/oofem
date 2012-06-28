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

#include "concretedpm.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "structuralms.h"
#include "gausspnt.h"
#include "intarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "structuralmaterial.h"
#include "isolinearelasticmaterial.h"
#include "structuralcrosssection.h"

namespace oofem {
ConcreteDPMStatus :: ConcreteDPMStatus(int n, Domain *d, GaussPoint *gp) :
    StructuralMaterialStatus(n, d, gp),
    plasticStrain( gp->giveMaterialMode() ),
    tempPlasticStrain( gp->giveMaterialMode() )
{
    kappaP = tempKappaP = 0.;
    kappaD = tempKappaD = 0.;
    damage = tempDamage = 0.;
    equivStrain = tempEquivStrain = 0.;
    deltaLambda = 0.;
    state_flag = temp_state_flag = ConcreteDPMStatus :: ConcreteDPM_Elastic;
#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    epsloc = -1.; // negative value means that it has not been set yet
#endif
}

ConcreteDPMStatus :: ~ConcreteDPMStatus()
{}

void
ConcreteDPMStatus :: initTempStatus()
{
    // Call the function of the parent class to initialize the variables defined there.
    StructuralMaterialStatus :: initTempStatus();
    tempPlasticStrain = plasticStrain;
    tempKappaP = kappaP;
    tempKappaD = kappaD;
    tempDamage = damage;
    tempEquivStrain = equivStrain;
    temp_state_flag = state_flag;
}

void
ConcreteDPMStatus :: updateYourself(TimeStep *atTime)
{
#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    // check whether the second-order work is negative
    // (must be done !!!before!!! the update of strain and stress)
    if ( epsloc < 0. ) {
        FloatArray strainIncrement = tempStrainVector;
        strainIncrement.subtract(strainVector);
        FloatArray stressIncrement = tempStressVector;
        stressIncrement.subtract(stressVector);
        int n = strainIncrement.giveSize();
        double work = strainIncrement.dotProduct(stressIncrement, n);
        //printf(" work : %g\n", work);
        if ( work < 0. ) {
            double E = gp->giveMaterial()->give('E', gp);
            double ft = gp->giveMaterial()->give(ft_strength, gp);
            epsloc = kappaD + damage * ft / E;
        }
    }

#endif

    // Call corresponding function of the parent class to update
    // variables defined there.
    StructuralMaterialStatus :: updateYourself(atTime);

    // update variables defined in ConcreteDPMStatus
    plasticStrain = tempPlasticStrain;
    kappaP = tempKappaP;
    kappaD = tempKappaD;
    damage = tempDamage;
    equivStrain = tempEquivStrain;
    state_flag = temp_state_flag;
}

void
ConcreteDPMStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    // Call corresponding function of the parent class to print
    // variables defined there.
    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "\tstatus { ");

    // print status flag
    switch ( state_flag ) {
    case ConcreteDPMStatus :: ConcreteDPM_Elastic:
        fprintf(file, "statusflag 0 (Elastic),");
        break;
    case ConcreteDPMStatus :: ConcreteDPM_Unloading:
        fprintf(file, "statusflag 1 (Unloading),");
        break;
    case ConcreteDPMStatus :: ConcreteDPM_Plastic:
        fprintf(file, "statusflag 2 (Plastic),");
        break;
    case ConcreteDPMStatus :: ConcreteDPM_Damage:
        fprintf(file, "statusflag 3 (Damage),");
        break;
    case ConcreteDPMStatus :: ConcreteDPM_PlasticDamage:
        fprintf(file, "statusflag 4 (PlasticDamage),");
        break;
    case ConcreteDPMStatus :: ConcreteDPM_VertexCompression:
        fprintf(file, "statusflag 5 (VertexCompression),");
        break;
    case ConcreteDPMStatus :: ConcreteDPM_VertexTension:
        fprintf(file, "statusflag 6 (VertexTension),");
        break;
    case ConcreteDPMStatus :: ConcreteDPM_VertexCompressionDamage:
        fprintf(file, "statusflag 7 (VertexCompressionDamage),");
        break;
    case ConcreteDPMStatus :: ConcreteDPM_VertexTensionDamage:
        fprintf(file, "statusflag 8 (VertexTensionDamage),");
        break;
    }

    /*
     * // print plastic strain vector
     * StrainVector plasticStrainVector( gp->giveMaterialMode() ) ;
     * giveFullPlasticStrainVector(plasticStrainVector) ;
     *
     * fprintf (file,"plastic strains ") ;
     * int n = plasticStrainVector.giveSize() ;
     * for (int i=1 ; i<=n ; i++)
     * fprintf (file," % .4e", plasticStrainVector.at(i)) ;
     */

    // print hardening/softening parameters and damage

    //fprintf (file," kappaP % .4e", kappaP ) ;
    //fprintf (file,", kappaD % .4e", kappaD ) ;
    fprintf(file, ", kappa % .4e % .4e", kappaP, kappaD);
    fprintf(file, ", damage % .4e", damage);

    // end of record
    fprintf(file, "}\n");
}

contextIOResultType
ConcreteDPMStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data

    if ( ( iores = plasticStrain.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->write(& kappaP, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& kappaD, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& equivStrain, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& state_flag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& deltaEquivStrain, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


contextIOResultType
ConcreteDPMStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( ( iores = plasticStrain.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->read(& kappaP, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& kappaD, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& equivStrain, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& state_flag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& deltaEquivStrain, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

int
ConcreteDPMStatus :: setIPValue(const FloatArray value, InternalStateType type)
{
    if ( type == IST_DamageScalar ) {
        damage = value.at(1);
        return 1;
    }

    if ( type == IST_CumPlasticStrain ) {
        kappaP = value.at(1);
        return 1;
    }

    if ( type == IST_CumPlasticStrain_2 ) {
        kappaD = value.at(1);
        return 1;
    }

    return 0;
}

/*
 * void
 * ConcreteDPMStatus::setStatusVariable (int varID, double value)
 * {
 * switch (varID){
 * case 1: damage = value; break;
 * case 2: kappaP = value; break;
 * case 3: kappaD = equivStrain = value; break;
 * case 4: stressVector.at(1) = value; break;
 * case 5: stressVector.at(2) = value; break;
 * case 6: stressVector.at(3) = value; break;
 * case 7: stressVector.at(4) = value; break;
 * case 8: stressVector.at(5) = value; break;
 * case 9: stressVector.at(6) = value; break;
 * default: printf("Warning: unknown variable ID %d in ConcreteDPMStatus::setStatusVariable\n",varID);
 * }
 * }
 */

void
ConcreteDPMStatus :: restoreConsistency()
{
    ConcreteDPM *mat = ( ConcreteDPM * ) gp->giveElement()->giveMaterial();

    // compute kappaD from damage
    kappaD = mat->computeInverseDamage(damage, gp);
    equivStrain = kappaD;

    // compute plastic strain
    // such that the given stress is obtained at zero total strain and given damage
    if ( damage < 1. ) {
        StressVector effectiveStress(stressVector, _3dMat);
        effectiveStress.times( -1. / ( 1. - damage ) );
        double E = mat->give('E', gp);
        double nu = mat->give('n', gp);
        effectiveStress.applyElasticCompliance(plasticStrain, E, nu);
    }
}


//   ******************************************************
//   *** CLASS CONCRETE DAMAGE-PLASTIC MATERIAL MODEL   ***
//   ******************************************************

//#define DPM_ITERATION_LIMIT 1.e-8

ConcreteDPM :: ConcreteDPM(int n, Domain *d) :
    StructuralMaterial(n, d),
    effectiveStress(_Unknown)
{
    tempKappaP = kappaP = 0.;
    tempKappaD = kappaD = 0.;
    tempDamage = damage = 0.;
    yieldTol = 0.;
    newtonIter = 0;
    matMode = _Unknown;
    helem = 0.;
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    href = 0.;
#endif
}

ConcreteDPM :: ~ConcreteDPM()
{
    delete linearElasticMaterial;
}

IRResultType
ConcreteDPM :: initializeFrom(InputRecord *ir)
{
    // Required by IR_GIVE_FIELD macro
    const char *__proc = "initializeFrom";
    IRResultType result;

    // call the corresponding service for the linear elastic material
    StructuralMaterial :: initializeFrom(ir);

    linearElasticMaterial->initializeFrom(ir);

    double value;
    // elastic parameters
    IR_GIVE_FIELD(ir, eM, IFT_IsotropicLinearElasticMaterial_e, "e");
    IR_GIVE_FIELD(ir, nu, IFT_IsotropicLinearElasticMaterial_n, "n");
    propertyDictionary->add('E', eM);
    propertyDictionary->add('n', nu);

    IR_GIVE_FIELD(ir, value, IFT_IsotropicLinearElasticMaterial_talpha, "talpha");
    propertyDictionary->add(tAlpha, value);

    gM = eM / ( 2. * ( 1. + nu ) );
    kM = eM / ( 3. * ( 1. - 2. * nu ) );

    // instanciate the variables of the plasticity model
    IR_GIVE_FIELD(ir, fc, IFT_ConcreteDPM_fc, "fc");
    IR_GIVE_FIELD(ir, ft, IFT_ConcreteDPM_ft, "ft");
    propertyDictionary->add(ft_strength, ft);
    propertyDictionary->add(fc_strength, fc);

    // damage parameters - only exponential softening
    // [in ef variable the wf (crack opening) is stored]
    if ( ir->hasField(IFT_ConcreteDPM_ef, "wf") ) {
        IR_GIVE_FIELD(ir, ef, IFT_ConcreteDPM_ef, "wf");
        // fracture energy
    } else {
        double Gf;
        IR_GIVE_FIELD(ir, Gf, IFT_ConcreteDPM_gf, "gf");
        ef = Gf / ft; // convert fracture energy to crack opening
    }

    // default parameters
    ecc = 0.525;
    IR_GIVE_OPTIONAL_FIELD(ir, ecc, IFT_ConcreteDPM_ecc, "ecc");
    yieldHardInitial = 0.1;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldHardInitial, IFT_ConcreteDPM_kinit, "kinit");
    AHard = 8.e-2;
    IR_GIVE_OPTIONAL_FIELD(ir, AHard, IFT_ConcreteDPM_ahard, "ahard");
    BHard = 3.e-3;
    IR_GIVE_OPTIONAL_FIELD(ir, BHard, IFT_ConcreteDPM_bhard, "bhard");
    CHard = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, CHard, IFT_ConcreteDPM_chard, "chard");
    DHard = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, DHard, IFT_ConcreteDPM_dhard, "dhard");
    ASoft = 15.;
    IR_GIVE_OPTIONAL_FIELD(ir, ASoft, IFT_ConcreteDPM_asoft, "asoft");
    helem = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, helem, IFT_ConcreteDPM_helem, "helem");
#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    href = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, href, IFT_ConcreteDPM_href, "href");
#endif

    //Compute m
    m = 3. * ( pow(fc, 2.) - pow(ft, 2.) ) / ( fc * ft ) * ecc / ( ecc + 1. );
    //Compute default value of dilationConst
    dilationConst = -0.85;
    IR_GIVE_OPTIONAL_FIELD(ir, dilationConst, IFT_ConcreteDPM_dilation, "dilation");

    yieldTol = 1.e-10;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, IFT_ConcreteDPM_yieldtol, "yieldtol");
    newtonIter = 100;
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, IFT_ConcreteDPM_newtoniter, "newtoniter");

    return IRRT_OK;
}

int
ConcreteDPM :: hasMaterialModeCapability(MaterialMode mMode)
{
    if ( ( mMode == _3dMat ) ||
        ( mMode == _PlaneStrain ) ) {
        return 1;
    } else {
        return 0;
    }
}

void
ConcreteDPM :: giveRealStressVector(FloatArray &answer,
                                    MatResponseForm form,
                                    GaussPoint *gp,
                                    const FloatArray &strainVector,
                                    TimeStep *atTime)
{
    FloatArray reducedTotalStrainVector;
    if ( matMode == _Unknown ) {
        matMode = gp->giveMaterialMode();
    }

    if ( effectiveStress.giveStressStrainMode() == _Unknown ) {
        effectiveStress.letStressStrainModeBe(matMode);
    }

    ConcreteDPMStatus *status = giveStatus(gp);

    // only for debugging
    //if (gp->giveElement()->giveNumber()==146){
    //  double a = 1.;
    //}

    // Initialize temp variables for this gauss point
    status->initTempStatus();

    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    // subtract stress-independent part of strain
    // (due to temperature changes, shrinkage, etc.)
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, strainVector, atTime, VM_Total);
    StrainVector strain( reducedTotalStrainVector, gp->giveMaterialMode() );

    // perform plasticity return
    performPlasticityReturn(gp, strain);

    // compute damage
    tempDamage = computeDamage(strain, gp, atTime);

    // compute elastic strains and effective stress
    StrainVector elasticStrain = strain;
    StrainVector tempPlasticStrain(matMode);
    status->giveTempPlasticStrain(tempPlasticStrain);
    elasticStrain.subtract(tempPlasticStrain);
    elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);

    // compute the nominal stress
    StressVector stress(matMode);
    stress = effectiveStress;

    stress.times(1. - tempDamage);

    status->letTempKappaDBe(tempKappaD);
    status->letTempDamageBe(tempDamage);

    status->letTempStrainVectorBe(strainVector);
    status->letTempStressVectorBe(stress);

    assignStateFlag(gp);

    if ( form == ReducedForm ) {
        answer = stress;
    } else {
        crossSection->giveFullCharacteristicVector(answer, gp, stress);
    }
}


double ConcreteDPM :: computeDamage(const StrainVector &strain, GaussPoint *gp, TimeStep *atTime)
{
    ConcreteDPMStatus *status = giveStatus(gp);
    double equivStrain;
    computeEquivalentStrain(equivStrain, strain, gp, atTime);
    double f = equivStrain - status->giveKappaD();
    if ( f <= 0.0 ) {
        // damage does not grow
        tempKappaD = status->giveKappaD();
        tempDamage = status->giveDamage();
    } else {
        // damage grows
        tempKappaD = equivStrain;
#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
        if ( ( href <= 0. ) || ( status->giveEpsLoc() >= 0. ) ) {
	  // evaluate and store the effective element size (if not known yet)
            this->initDamaged(tempKappaD, strain, gp);
        }

#else
        this->initDamaged(tempKappaD, strain, gp);
#endif
        tempDamage = computeDamageParam(tempKappaD, gp);
    }

    return tempDamage;
}

void
ConcreteDPM :: computeEquivalentStrain(double &tempEquivStrain, const StrainVector &strain, GaussPoint *gp, TimeStep *atTime)
{
    //The equivalent  strain is based on the volumetric plastic strain
    ConcreteDPMStatus *status = giveStatus(gp);
    MaterialMode matMode = gp->giveMaterialMode();
    tempKappaP = status->giveTempKappaP();
    kappaP = status->giveKappaP();
    double equivStrain = status->giveEquivStrain();
    double deltaEquivStrain = 0.;
    if ( tempKappaP <= 1.0 || tempKappaP == kappaP ) {
        tempEquivStrain = equivStrain;
        return;
    } else if ( tempKappaP > 1.0 && tempKappaP != kappaP ) {
        StrainVector plasticStrain(matMode);
        StrainVector tempPlasticStrain(matMode);
        status->giveTempPlasticStrain(tempPlasticStrain);
        status->givePlasticStrain(plasticStrain);
        double volumetricPlasticStrain = plasticStrain(0) + plasticStrain(1) +
                                         plasticStrain(2);
        double tempVolumetricPlasticStrain = tempPlasticStrain(0) +
                                             tempPlasticStrain(1) + tempPlasticStrain(2);
        if ( kappaP < 1.0 ) {
            //compute volumetric plastic strain at peak
            double peakVolumetricPlasticStrain = ( 1. - kappaP ) / ( tempKappaP - kappaP ) *
                                                 ( tempVolumetricPlasticStrain - volumetricPlasticStrain ) +
                                                 volumetricPlasticStrain;
            if ( peakVolumetricPlasticStrain < 0. ) {
                peakVolumetricPlasticStrain = 0.;
            }

            deltaEquivStrain =
                tempVolumetricPlasticStrain - peakVolumetricPlasticStrain;
            tempEquivStrain =
                deltaEquivStrain / computeDuctilityMeasureDamage(strain, gp);
            if ( tempEquivStrain < 0. ) {
                tempEquivStrain = 0.;
            }
        } else {
            deltaEquivStrain =  ( tempVolumetricPlasticStrain - volumetricPlasticStrain );
            if ( deltaEquivStrain < 0. ) {
                deltaEquivStrain = 0.;
            }

            tempEquivStrain = equivStrain +
                              deltaEquivStrain / computeDuctilityMeasureDamage(strain, gp);
        }
    }

    status->letTempEquivStrainBe(tempEquivStrain);
    status->letDeltaEquivStrainBe(deltaEquivStrain);
}

double
ConcreteDPM :: computeInverseDamage(double dam, GaussPoint *gp)
{
    ConcreteDPMStatus *status = giveStatus(gp);
    double le = status->giveLe();
    if ( le == 0. ) {
        if ( helem > 0. ) {
            le = helem;
        } else {
            le = gp->giveElement()->computeMeanSize();
        }

        status->setLe(helem);
    }

    double answer = -( this->ef / le ) * log(1. - dam) - dam * ft / eM;
    return answer;
}

#define DPM_DAMAGE_TOLERANCE 1.e-8

double
ConcreteDPM :: computeDamageParam(double kappa, GaussPoint *gp)
{
    double omega = 0.;
    if ( kappa > 0. ) {
        // iteration to achieve mesh objectivity
        // this is only valid for tension
        ConcreteDPMStatus *status = giveStatus(gp);
        int nite = 0;
        double R, Lhs, aux, aux1;

	double h = status->giveLe(); // effective element size

#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
	// the procedure is different before and after localization
        double epsloc = status->giveEpsLoc();
        if ( ( href <= 0. ) || ( epsloc < 0. ) ) { // before localization, the reference element size href is used (but only if it is really specified by the user)
	  if ( href > 0. ) {
	    h = href; // reference element size
	  }
#endif
	  // standard damage evaluation procedure
	  aux1 = ( ft / eM ) * h / ef;
	  if ( aux1 > 1 ) {
            printf("computeDamageParam: ft=%g, E=%g, wf=%g, hmax=E*wf/ft=%g, h=%g\n", ft, eM, ef, eM * ef / ft, h);
            _error("computeDamageParam: element too large");
	  }
	  do {
            nite++;
            aux = exp(-h * ( omega * ft / eM + kappa ) / ef);
            R = 1. - omega - aux;
            Lhs = -1. + aux1 * aux;
            omega -= R / Lhs;
            if ( nite > 40 ) {
	      _error("computeDamageParam: algorithm not converging");
            }
	  }
	  while ( fabs(R) >= DPM_DAMAGE_TOLERANCE );

#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
	} else   { // after localization, more complicated formula
	  if ( helem > 0. ) {
            h = helem; // predefined element size
	  } else {
            h = status->giveLe(); // effective element size
	  }

	  aux1 = ( ft / eM ) * h / ef;
	  do {
            nite++;
            aux = exp(-( ( href - h ) * epsloc + h * ( omega * ft / eM + kappa ) ) / ef);
            R = 1. - omega - aux;
            Lhs = -1. + aux1 * aux;
            omega -= R / Lhs;
            if ( nite > 40 ) {
	      printf("computeDamageParam: algorithm not converging (part 2)");
            }
	  } 
	  while ( fabs(R) >= DPM_DAMAGE_TOLERANCE );
	}
#endif
        if ( ( omega > 1.0 ) || ( omega < 0.0 ) ) {
	  _error2("computeDamageParam: internal error, omega = %g", omega);
        }
    }

    return omega;
}


void
ConcreteDPM :: initDamaged(double kappaD,
                           const StrainVector &strain,
                           GaussPoint *gp)
{
    if ( kappaD <= 0. ) {
        return;
    }

    int i, indx = 1;
    double le;
    FloatArray principalStrains, crackPlaneNormal(3);
    FloatMatrix principalDir(3, 3);
    ConcreteDPMStatus *status = giveStatus(gp);

    if ( helem > 0. ) {
      status->setLe(helem);
    } else if ( status->giveDamage() == 0. ) {
      strain.computePrincipalValDir(principalStrains, principalDir);
      // find index of max positive principal strain
      for ( i = 2; i <= 3; i++ ) {
        if ( principalStrains.at(i) > principalStrains.at(indx) ) {
	  indx = i;
        }
      }

      for ( i = 1; i <= 3; i++ ) {
        crackPlaneNormal.at(i) = principalDir.at(i, indx);
      }

    // Warning: This would give the element size divided by the number of Gauss points
    // le = gp->giveElement()->giveCharacteristicLenght (gp, crackPlaneNormal);

    // this gives the projected element size
    le = gp->giveElement()->giveLenghtInDir(crackPlaneNormal);
    if ( le == 0. ) {
      le = gp->giveElement()->computeMeanSize();
    }
    // store le in the corresponding status
    status->setLe(le);
    } else if ( status->giveLe() == 0. ) {
      // this happens if the status is initialized from a file
      // with nonzero damage
      // le determined as square root of element area or cube root of el. volume
      le = gp->giveElement()->computeMeanSize();
      status->setLe(le);
    }
}

double
ConcreteDPM :: computeDuctilityMeasureDamage(const StrainVector &strain, GaussPoint *gp)
{
    ConcreteDPMStatus *status = giveStatus(gp);
    StrainVector plasticStrain(matMode);
    StrainVector tempPlasticStrain(matMode);
    status->giveTempPlasticStrain(tempPlasticStrain);
    status->givePlasticStrain(plasticStrain);
    tempPlasticStrain.subtract(plasticStrain);
    StrainVector principalStrain(matMode);
    double ductilityMeasure;
    double volStrain = tempPlasticStrain(0) + tempPlasticStrain(1) +
                       tempPlasticStrain(2);
    //compute sum of negative principal strains
    tempPlasticStrain.computePrincipalValues(principalStrain);
    double negativeVolStrain = 0.;
    for ( int i = 0; i < principalStrain.giveSize(); i++ ) {
        if ( principalStrain(i) < 0. ) {
            negativeVolStrain += principalStrain(i);
        }
    }

    //compute ductility measure
    double Rs = -negativeVolStrain / volStrain;
    if ( Rs < 1.0 ) {
        ductilityMeasure = 1. + ASoft *pow(Rs, 2.);
    } else {
        ductilityMeasure = 1. - 3. * ASoft + 4. *ASoft *sqrt(Rs);
    }

    return ductilityMeasure;
}

void
ConcreteDPM :: performPlasticityReturn(GaussPoint *gp,
                                       StrainVector &strain)
{
    ConcreteDPMStatus *status = ( ConcreteDPMStatus * ) this->giveStatus(gp);
    if ( matMode == _Unknown ) {
        matMode = gp->giveMaterialMode();
    }

    if ( effectiveStress.giveStressStrainMode() == _Unknown ) {
        effectiveStress.letStressStrainModeBe(matMode);
    }

    status->initTempStatus();

    //get temp plastic strain and tempKappa
    StrainVector tempPlasticStrain(matMode);
    status->giveTempPlasticStrain(tempPlasticStrain);
    tempKappaP = status->giveTempKappaP();

    // compute elastic strains and trial stress
    StrainVector elasticStrain = strain;
    elasticStrain.subtract(tempPlasticStrain);
    elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);

    //Compute trial coordinates
    computeTrialCoordinates(effectiveStress);

    double yieldValue = computeYieldValue(sig, rho, thetaTrial, tempKappaP);
    // choose correct stress return and update state flag
    if ( yieldValue  > yieldTol ) {
        double apexStress = 0.;
        checkForVertexCase(apexStress, sig, tempKappaP);

        //Make the appropriate return
        if ( vertexType == VT_Tension || vertexType == VT_Compression ) {
            performVertexReturn(effectiveStress, apexStress, gp);
            if ( vertexType == VT_Regular ) {
                //This was no real vertex case
                //get the original tempKappaP and stress
                tempKappaP = status->giveTempKappaP();
                elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);
            }
        }

        if ( vertexType == VT_Regular ) {
            performRegularReturn(effectiveStress, gp);
        }
    }

    // update temp kappaP
    status->letTempKappaPBe(tempKappaP);
    // compute the plastic strains
    effectiveStress.applyElasticCompliance(elasticStrain, eM, nu);
    tempPlasticStrain = strain;
    tempPlasticStrain.subtract(elasticStrain);
    status->letTempPlasticStrainBe(tempPlasticStrain);
}


bool
ConcreteDPM :: checkForVertexCase(double answer,
                                  const double sig,
                                  const double tempKappa)
{
    //Compute sigZero and compare with actual sig
    int l = 0;
    double FZero = 1.;
    double dFZeroDSigZero = 0.;
    double sigZero;
    const double yieldHardOne = computeHardeningOne(tempKappa);
    double apexCompression = 0.;

    //compressive apex
    if ( ( tempKappa < 1. ) && ( sig < 0. ) ) {
        sigZero = -15. * fc;
        while ( fabs(FZero) > yieldTol && l <= newtonIter ) {
            l++;
            FZero = pow( ( 1. - yieldHardOne ), 2. ) * pow( ( sigZero / fc ), 4. ) +
                    pow(yieldHardOne, 2.) * m * ( sigZero / fc ) - pow(yieldHardOne, 2.);

            dFZeroDSigZero = pow( ( 1. - yieldHardOne ), 2. ) * 4. * pow( ( sigZero / fc ), 3. ) / fc +
                             pow(yieldHardOne, 2.) * m / fc;

            sigZero = sigZero - FZero / dFZeroDSigZero;
        }

        if ( l < 15 && sigZero < 0. ) {
            apexCompression = sigZero;
        } else {
            apexCompression = -15. * fc;
        }
    }

    if ( ( sig > 0. && tempKappa < 1. ) || ( sig > fc / m && tempKappa >= 1. ) ) {
        vertexType = VT_Tension;
        answer = 0.;
    } else if ( sig < apexCompression ) {
        vertexType = VT_Compression;
        answer = apexCompression;
    } else {
        vertexType = VT_Regular;
    }

    return false;
}



void
ConcreteDPM :: performRegularReturn(StressVector &effectiveStress,
                                    GaussPoint *gp)
{
    //Define status
    ConcreteDPMStatus *status = ( ConcreteDPMStatus * ) this->giveStatus(gp);
    //Variables
    deltaLambda = 0.;
    double deltaLambdaIncrement = 0.;
    double yieldValue = 1.;
    double residualNorm = 1.;
    double rhoTest = 0.;
    double deltaLambdaIncrementNew = 0.;
    double tempKappaPTest = 0.;
    double deltaKappaP = 0.;
    int iterationCount = 0;
    int negativeRhoFlag = 0;
    StrainVector elasticStrain(matMode);
    FloatArray dGDInv(2);
    FloatArray dFDInv(2);
    FloatArray residual(3);
    FloatArray derivativesOfYieldSurface(3);
    FloatArray flowRules(3);
    FloatArray answerIncrement(3);
    //  double residual2;
    double dKappaDDeltaLambda;
    double dFDKappa = 0.;
    FloatMatrix aMatrix;
    FloatArray helpVectorA(3);
    FloatArray helpVectorB(3);
    StressVector helpIncrement(matMode);
    StrainVector plasticStrainIncrement(matMode);
    StrainVector dGDStressDeviatoric(matMode);
    StressVector deviatoricStress(matMode);

    //compute the principal directions of the stress
    FloatArray help;
    FloatMatrix stressPrincipalDir;
    effectiveStress.computePrincipalValDir(help, stressPrincipalDir);

    //compute invariants from stress state
    effectiveStress.computeDeviatoricVolumetricSplit(deviatoricStress, sig);
    rho = deviatoricStress.computeSecondCoordinate();

    const StressVector initialStress = effectiveStress;

    const double volumetricPlasticStrain =
        3. * status->giveVolumetricPlasticStrain();

    const double deviatoricPlasticStrainNorm =
        status->giveDeviatoricPlasticStrainNorm();

    double tempVolumetricPlasticStrain = volumetricPlasticStrain;
    double tempDeviatoricPlasticStrainNorm = deviatoricPlasticStrainNorm;

    const double kappaP = status->giveKappaP();
    tempKappaP = kappaP;

    yieldValue = computeYieldValue(sig, rho, thetaTrial, tempKappaP);
    residualNorm = fabs(yieldValue);

    while ( ( residualNorm ) > yieldTol ) {
        if ( ++iterationCount == newtonIter ) {
            _error("Closest point projection did not converge.\n");
        }

        //compute the stress, yield value and residuals
        computeDGDInv(dGDInv, sig, rho, tempKappaP);
        computeDFDInv(dFDInv, sig, rho, tempKappaP);
        dFDKappa = computeDFDKappa(sig, rho, tempKappaP);
        dKappaDDeltaLambda = computeDKappaDDeltaLambda(sig, rho, tempKappaP);
        yieldValue = computeYieldValue(sig, rho, thetaTrial, tempKappaP);
        //The residual volumetric plastic strain is defined as eps1+eps2+eps3
        residual(0) = volumetricPlasticStrain - tempVolumetricPlasticStrain +
                      deltaLambda *dGDInv(0);
        residual(1) = deviatoricPlasticStrainNorm -
                      tempDeviatoricPlasticStrainNorm + deltaLambda *dGDInv(1);
        residual(2) = kappaP - tempKappaP + deltaLambda * dKappaDDeltaLambda;

        // weighted norm
        residualNorm = pow(residual(0), 2.) + pow(residual(1), 2.) +
                       pow(residual(2), 2.) + pow(yieldValue, 2.);
        residualNorm = sqrt(residualNorm);
        //    printf("\n residualNorm= %e\n", residualNorm);
        if ( residualNorm > yieldTol ) {
            computeAMatrix(aMatrix, sig, rho, tempKappaP);

            // assemble the derivatives of the yield surface
            for ( int i = 0; i < 2; i++ ) {
                derivativesOfYieldSurface(i) = dFDInv(i);
            }

            derivativesOfYieldSurface(2) = dFDKappa;
            //assemble flow rules
            for ( int i = 0; i < 2; i++ ) {
                flowRules(i) = dGDInv(i);
            }

            flowRules(2) = dKappaDDeltaLambda;

            //compute the deltaLambdaIncrement
            deltaLambdaIncrement = yieldValue;
            helpVectorA.beProductOf(aMatrix, residual);
            for ( int i = 0; i < 3; i++ ) {
                deltaLambdaIncrement -= derivativesOfYieldSurface(i) * helpVectorA(i);
            }

            helpVectorB.beProductOf(aMatrix, flowRules);
            double denominator = 0.;
            for ( int i = 0; i < 3; i++ ) {
                denominator += derivativesOfYieldSurface(i) * helpVectorB(i);
            }

            deltaLambdaIncrement /= denominator;

            //compute increment of answer
            FloatArray helpVectorC;
            helpVectorC = residual;
            helpVectorC.add(deltaLambdaIncrement, flowRules);
            answerIncrement.beProductOf(aMatrix, helpVectorC);
            answerIncrement.negated();
            rhoTest = rho + answerIncrement(1);

            //Special case, if rho changes sign
            if ( rhoTest < 0. && negativeRhoFlag == 0 ) {
                //Determine deltaLambdaIncrement, so that rho is equal to zero
                answerIncrement(1) = -rho;

                deltaLambdaIncrementNew =
                    ( -aMatrix(1, 0) * residual(0) - aMatrix(1, 1) * residual(1) -
                     aMatrix(1, 2) * residual(2) - answerIncrement(1) ) /
                    ( flowRules(0) * aMatrix(1, 0) + flowRules(1) * aMatrix(1, 1) +
                     flowRules(2) * aMatrix(1, 2) );

                // Special case, if deltaLambdaIncrement is equal to zero.
                if ( fabs(deltaLambdaIncrementNew) < yieldTol * 1.e3 ) {
                    negativeRhoFlag = 1;
                    deltaLambdaIncrementNew = deltaLambdaIncrement;
                }

                helpVectorC = residual;
                helpVectorC.add(deltaLambdaIncrementNew, flowRules);
                answerIncrement.beProductOf(aMatrix, helpVectorC);
                answerIncrement.negated();

                sig += answerIncrement(0);
                rho += answerIncrement(1);

                tempKappaPTest = tempKappaP;
                tempKappaP += answerIncrement(2);
                deltaLambda += deltaLambdaIncrementNew;
            } else {
                tempKappaPTest = tempKappaP;
                tempKappaP += answerIncrement(2);
                deltaKappaP = tempKappaP - status->giveKappaP();
                sig += answerIncrement(0);
                rho += answerIncrement(1);

                deltaLambda += deltaLambdaIncrement;
            }

            //Special case, if deltaKappaP is negative
            if ( ( tempKappaP - status->giveKappaP() ) < 0. ) {
                tempKappaP = tempKappaPTest;
            }

            tempVolumetricPlasticStrain -= answerIncrement(0) / ( kM );
            tempDeviatoricPlasticStrainNorm -= answerIncrement(1) / ( 2. * gM );
        }

        //if (gp->giveElement()->giveNumber() == 1301){
        //  printf("%g %g %g\n",sig,rho,tempKappaP);
        //}
    }

    status->letTempVolumetricPlasticStrainBe(tempVolumetricPlasticStrain / 3.);
    if ( deltaLambda < 0 ) {
        printf("deltaLambda = %e\n", deltaLambda);
        printf("plastic multiplier less than zero");
    }

    // printf("\nnumber of iterations = %i\n", iterationCount);
    status->letDeltaLambdaBe(deltaLambda);
    StressVector stressPrincipal(_3dMat);
    StressVector stressTemp(_3dMat);
    stressPrincipal.zero();

    stressPrincipal(0) = sig + sqrt(2. / 3.) * rho * cos(thetaTrial);
    stressPrincipal(1) = sig + sqrt(2. / 3.) * rho * cos(thetaTrial - 2. * M_PI / 3.);
    stressPrincipal(2) = sig + sqrt(2. / 3.) * rho * cos(thetaTrial + 2. * M_PI / 3.);

    transformStressVectorTo(stressTemp, stressPrincipalDir, stressPrincipal, 1);

    if ( matMode == _PlaneStrain ) {
        effectiveStress(0) = stressTemp(0);
        effectiveStress(1) = stressTemp(1);
        effectiveStress(2) = stressTemp(2);
        effectiveStress(3) = stressTemp(5);
    } else {
        effectiveStress = stressTemp;
    }
}

void
ConcreteDPM :: performVertexReturn(StressVector &effectiveStress,
                                   double apexStress,
                                   GaussPoint *gp)
{
    ConcreteDPMStatus *status = ( ConcreteDPMStatus * ) this->giveStatus(gp);
    StressVector deviatoricStressTrial(matMode);
    double sigTrial;
    StressVector stressTemp(matMode);
    StressVector deviatoricStress(matMode);
    double yieldValue = 0.;
    double yieldValueMid = 0.;
    double sig2 = 0.;
    double dSig;
    double sigMid;
    double sigAnswer;
    double ratioPotential;

    effectiveStress.computeDeviatoricVolumetricSplit(deviatoricStressTrial, sigTrial);
    const double rhoTrial = deviatoricStressTrial.computeSecondCoordinate();

    StrainVector deltaPlasticStrain(matMode);

    tempKappaP = status->giveTempKappaP();
    const double kappaInitial = tempKappaP;

    if ( vertexType == VT_Tension ) {
        sig2 = -0.1 * ft;
    } else if ( vertexType == VT_Compression ) {
        sig2 = apexStress;
    }

    tempKappaP =
        computeTempKappa(kappaInitial, sigTrial, rhoTrial, sigTrial);

    yieldValue =
        computeYieldValue(sigTrial, 0., 0., tempKappaP);

    tempKappaP =
        computeTempKappa(kappaInitial, sigTrial, rhoTrial, sig2);

    yieldValueMid =
        computeYieldValue(sig2, 0., 0., tempKappaP);

    if ( yieldValue * yieldValueMid >= 0. ) {
        vertexType = VT_Regular;
        return;
    }

    if ( yieldValue < 0.0 ) {
        dSig = sig2 - sigTrial;
        sigAnswer = sig2;
    } else {
        dSig = sigTrial - sig2;
        sigAnswer = sig2;
    }

    for ( int j = 0; j < 250; j++ ) {
        dSig = 0.5 * dSig;

        sigMid = sigAnswer + dSig;


        tempKappaP =
            computeTempKappa(kappaInitial, sigTrial, rhoTrial, sigMid);

        yieldValueMid = computeYieldValue(sigMid, 0., 0., tempKappaP);

        if ( yieldValueMid <= 0. ) {
            sigAnswer = sigMid;
        }

        if ( fabs(yieldValueMid) < yieldTol && yieldValueMid <= 0. ) {
            for ( int i = 0; i < 3; i++ ) {
                effectiveStress(i) = sigAnswer;
            }

            for ( int i = 3; i < effectiveStress.giveSize(); i++ ) {
                effectiveStress(i) = 0.;
            }

            ratioPotential =
                computeRatioPotential(sigAnswer, tempKappaP);

            double ratioTrial = rhoTrial / ( sigTrial - sigAnswer );

            if ( ( ( ( ratioPotential >= ratioTrial ) && vertexType == VT_Tension ) ) ||
                ( ( ratioPotential <= ratioTrial ) && vertexType == VT_Compression ) ) {
                return;
            } else {
                vertexType = VT_Regular;
                return;
            }
        }
    }

    vertexType = VT_Regular;
}


double
ConcreteDPM :: computeTempKappa(const double kappaInitial,
                                const double sigTrial,
                                const double rhoTrial,
                                const double sig)
{
    //This function is called, if stress state is in vertex case
    double equivalentDeltaPlasticStrain;
    FloatArray deltaPlasticStrainPrincipal(3);
    rho = 0.;
    equivalentDeltaPlasticStrain = sqrt( 1. / 9. * pow( ( sigTrial - sig ) / ( kM ), 2. ) +
                                        pow(rhoTrial / ( 2. * gM ), 2.) );

    double ductilityMeasure = computeDuctilityMeasure(sig, rho);

    return kappaInitial + equivalentDeltaPlasticStrain / ductilityMeasure;
}


double
ConcreteDPM :: computeYieldValue(const double sig,
                                 const double rho,
                                 const double theta,
                                 const double tempKappa) const
{
    //compute yieldHard
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double yieldHardTwo = computeHardeningOne(tempKappa);

    //  compute elliptic function r
    const double rFunction = ( 4. * ( 1. - pow(ecc, 2.) ) * pow(cos(theta), 2.) +
                              pow( ( 2. * ecc - 1. ), 2. ) ) /
                             ( 2. * ( 1. - pow(ecc, 2.) ) * cos(theta) +
                              ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - pow(ecc, 2.) ) * pow(cos(theta), 2.)
                                                       + 5. * pow(ecc, 2.) - 4. * ecc) );

    //compute help function Al
    const double Al = ( 1. - yieldHardOne ) * pow( ( sig / fc + rho / ( sqrt(6.) * fc ) ), 2. ) +
                      sqrt(3. / 2.) * rho / fc;

    //Compute yield equation
    return pow(Al, 2.) +
           pow(yieldHardOne, 2.) * m * ( sig / fc + rho * rFunction / ( sqrt(6.) * fc ) ) -
           pow(yieldHardTwo, 2.);
}


double
ConcreteDPM :: computeDFDKappa(const double sig,
                               const double rho,
                               const double tempKappa)
{
    const double theta = thetaTrial;
    //compute yieldHard and yieldSoft
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double yieldHardTwo = computeHardeningOne(tempKappa);
    // compute the derivative of the hardening and softening laws
    const double dYieldHardOneDKappa = computeHardeningOnePrime(tempKappa);
    const double dYieldHardTwoDKappa = computeHardeningOnePrime(tempKappa);

    //compute elliptic function r
    const double rFunction =
        ( 4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta) + ( 2. * ecc - 1. ) * ( 2. * ecc - 1. ) ) /
        ( 2 * ( 1. - ecc * ecc ) * cos(theta) + ( 2. * ecc - 1. ) *
         sqrt(4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta) + 5. * ecc * ecc - 4. * ecc) );

    //compute help functions Al, Bl
    const double Al = ( 1. - yieldHardOne ) * pow( ( sig / fc + rho / ( sqrt(6.) * fc ) ), 2. ) + sqrt(3. / 2.) * rho / fc;

    const double Bl = sig / fc + rho / ( fc * sqrt(6.) );

    const double dFDYieldHardOne = -2. *Al *pow(Bl, 2.)
                                   + 2. * yieldHardOne * m * ( sig / fc + rho * rFunction / ( sqrt(6.) * fc ) );

    const double dFDYieldHardTwo = -2. * yieldHardTwo;

    // compute dFDKappa
    double dFDKappa =  dFDYieldHardOne * dYieldHardOneDKappa +
                      dFDYieldHardTwo * dYieldHardTwoDKappa;
    /*
     * set dFDKappa to zero, if it becomes greater than zero.
     * dFDKappa can only be negative or zero in the converged state for
     * the case of hardenig and perfect plasticity. For trial stresses, however,
     * it might be negative, which may lead to convergency problems. Therefore,
     * it is set to zero in this cases.
     */
    if ( dFDKappa > 0. ) {
        dFDKappa = 0.;
    }

    return dFDKappa;
}

void
ConcreteDPM :: computeDFDInv(FloatArray &answer,
                             const double sig,
                             const double rho,
                             const double tempKappa) const
{
    const double theta = thetaTrial;

    //compute yieldHard
    const double yieldHardOne = computeHardeningOne(tempKappa);

    //compute elliptic function r
    const double rFunction = ( 4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta) + ( 2. * ecc - 1. ) * ( 2. * ecc - 1. ) ) /
                             ( 2. * ( 1. - ecc * ecc ) * cos(theta) + ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta)
                                                                                               + 5. * ecc * ecc - 4. * ecc) );

    //compute help functions AL, BL
    const double AL = ( 1. - yieldHardOne ) * pow( ( sig / fc + rho / ( sqrt(6.) * fc ) ), 2. ) + sqrt(3. / 2.) * rho / fc;
    const double BL = sig / fc + rho / ( fc * sqrt(6.) );

    //compute dfdsig
    const double dfdsig = 4. * ( 1. - yieldHardOne ) / fc * AL * BL + yieldHardOne * yieldHardOne * m / fc;
    //compute dfdrho
    const double dfdrho = AL / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * BL + 6. ) + rFunction * m * yieldHardOne * yieldHardOne / ( sqrt(6.) * fc );

    answer(0) = dfdsig;
    answer(1) = dfdrho;
}

double
ConcreteDPM :: computeDKappaDDeltaLambda(const double sig,
                                         const double rho,
                                         const double tempKappa)
{
    //Variables
    double equivalentDGDStress;
    FloatArray dGDInv(2);
    computeDGDInv(dGDInv, sig, rho, tempKappa);
    FloatArray dGDStressPrincipal(3);

    equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv(0), 2.) +
                               pow(dGDInv(1), 2.) );

    double ductilityMeasure = computeDuctilityMeasure(sig, rho);
    double dKappaDDeltaLambda = equivalentDGDStress / ductilityMeasure;
    return dKappaDDeltaLambda;
}


void
ConcreteDPM :: computeDDKappaDDeltaLambdaDInv(FloatArray &answer,
                                              const double sig,
                                              const double rho,
                                              const double tempKappa)
{
    //Variables
    double equivalentDGDStress;
    FloatArray dGDInv(2);
    FloatMatrix dDGDDInv(2, 2);
    FloatArray dGDStressPrincipal(3);
    FloatArray dEquivalentDGDStressDInv(2);
    FloatArray helpA(3);

    //Compute first and second derivative of plastic potential
    computeDGDInv(dGDInv, sig, rho, tempKappa);
    computeDDGDDInv(dDGDDInv, sig, rho, tempKappa);

    //Compute equivalentDGDStress
    equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv(0), 2.) +
                               pow(dGDInv(1), 2.) );

    //computeDuctilityMeasure
    double ductilityMeasure = computeDuctilityMeasure(sig, rho);

    //Compute dEquivalentDGDStressDInv
    dEquivalentDGDStressDInv(0) =
        ( 2. / 3. * dGDInv(0) * dDGDDInv(0, 0) + 2. * dGDInv(1) * dDGDDInv(1, 0) ) / ( 2. * equivalentDGDStress );
    dEquivalentDGDStressDInv(1) =
        ( 2. / 3. * dGDInv(0) * dDGDDInv(0, 1) + 2. * dGDInv(1) * dDGDDInv(1, 1) ) / ( 2. * equivalentDGDStress );

    answer.zero();

    // compute the derivative of
    FloatArray dDuctilityMeasureDInv(2);
    computeDDuctilityMeasureDInv(dDuctilityMeasureDInv, sig, rho, tempKappa);

    answer(0) = ( dEquivalentDGDStressDInv(0) * ductilityMeasure - equivalentDGDStress * dDuctilityMeasureDInv(0) ) / pow(ductilityMeasure, 2.);

    answer(1) = ( dEquivalentDGDStressDInv(1) * ductilityMeasure - equivalentDGDStress * dDuctilityMeasureDInv(1) ) / pow(ductilityMeasure, 2.);
}

double
ConcreteDPM :: computeDDKappaDDeltaLambdaDKappa(const double sig,
                                                const double rho,
                                                const double tempKappa)
{
    //Variables
    double equivalentDGDStress;
    FloatArray dGDInv(2);
    FloatArray dDGDInvDKappa(2);
    FloatArray dGDStressPrincipal(3);
    FloatArray helpA(3);
    double dEquivalentDGDStressDKappa;

    //Compute first and second derivative of plastic potential
    computeDGDInv(dGDInv, sig, rho, tempKappa);
    computeDDGDInvDKappa(dDGDInvDKappa, sig, rho, tempKappa);

    equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv(0), 2.) +
                               pow(dGDInv(1), 2.) );

    //computeDuctilityMeasure
    double ductilityMeasure = computeDuctilityMeasure(sig, rho);

    //Compute dEquivalentDGDStressDKappa
    dEquivalentDGDStressDKappa =
        ( 2. / 3. * dGDInv(0) * dDGDInvDKappa(0) + 2. * dGDInv(1) * dDGDInvDKappa(1) ) / ( 2. * equivalentDGDStress );

    // compute the derivative of
    double dDuctilityMeasureDKappa = 0.;

    double dDKappaDDeltaLambdaDKappa =
        ( dEquivalentDGDStressDKappa * ductilityMeasure -
          equivalentDGDStress * dDuctilityMeasureDKappa ) / pow(ductilityMeasure, 2.);

    return dDKappaDDeltaLambdaDKappa;
}

double
ConcreteDPM :: computeDuctilityMeasure(const double sig,
                                       const double rho)
{
    double thetaConst = pow(2. * cos(thetaTrial), 2.);
    double ductilityMeasure;
    double x = -( sig + fc / 3 ) / fc;
    if ( x < 0. ) {
        /*Introduce exponential help function which results in a smooth
         * transition. */
        double fZero = BHard;
        double fPrimeZero = -( BHard - AHard ) / ( CHard );
        double CHelp = DHard;
        double AHelp = fZero - CHelp;
        double BHelp = ( fZero - CHelp ) / fPrimeZero;
        ductilityMeasure = ( AHelp * exp(x / BHelp) + CHelp ) / thetaConst;
    } else {
        ductilityMeasure = ( AHard + ( BHard - AHard ) * exp( -x / ( CHard ) ) ) / thetaConst;
    }

    if ( ductilityMeasure <= 0. ) {
        printf("ductilityMeasure is zero or negative");
    }

    return ductilityMeasure;
}

void
ConcreteDPM :: computeDDuctilityMeasureDInv(FloatArray &answer,
                                            const double sig,
                                            const double rho,
                                            const double tempKappa)
{
    double thetaConst = pow(2. * cos(thetaTrial), 2.);
    double x = ( -( sig + fc / 3 ) ) / fc;
    if ( x < 0. ) {
        double dXDSig = -1. / fc;
        /* Introduce exponential help function which results in a
         * smooth transition. */
        double fZero = BHard;
        double fPrimeZero =  -( BHard - AHard ) / ( CHard );
        double CHelp = DHard;
        double AHelp = fZero - CHelp;
        double BHelp = ( fZero - CHelp ) / fPrimeZero;
        double dDuctilityMeasureDX = AHelp / BHelp *exp(x / BHelp) / thetaConst;
        answer(0) = dDuctilityMeasureDX * dXDSig;
        answer(1) = 0.;
    } else {
        double dXDSig = -1. / fc;
        double dDuctilityMeasureDX = -( BHard - AHard ) / ( CHard ) / thetaConst *exp( -x / ( CHard ) );
        answer(0) = dDuctilityMeasureDX * dXDSig;
        answer(1) = 0.;
    }
}


void
ConcreteDPM :: computeDGDInv(FloatArray &answer,
                             const double sig,
                             const double rho,
                             const double tempKappa)
{
    //Compute dilation parameter
    const double R = ( sig - ft / 3. ) / fc;
    const double mQTension = ( 3. * ft / fc + m / 2. );
    const double mQCompression =
        ( 1. + 2. * dilationConst ) / ( dilationConst - 1. ) * ( 3. + m / 2. );
    const double AConst = -( ft + fc ) / ( 3. * fc ) / log(mQCompression / mQTension);
    const double mQ = mQTension * exp(R / AConst);

    //compute yieldHard and yieldSoft
    const double yieldHardOne = computeHardeningOne(tempKappa);

    const double Bl = sig / fc + rho / ( fc * sqrt(6.) );

    const double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fc;

    const double dgdsig = 4. * ( 1. - yieldHardOne ) / fc * Al * Bl + yieldHardOne * yieldHardOne * mQ / fc;

    const double dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +
                          m *pow(yieldHardOne, 2.) / ( sqrt(6.) * fc );

    answer(0) = dgdsig;
    answer(1) = dgdrho;
}


double
ConcreteDPM :: computeRatioPotential(const double sig,
                                     const double tempKappa)
{
    //Compute dilation parameter
    const double R = ( sig - ft / 3. ) / fc;
    const double mQTension = ( 3. * ft / fc + m / 2. );
    const double mQCompression =
        ( 1. + 2. * dilationConst ) / ( dilationConst - 1. ) * ( 3. + m / 2. );
    const double AConst = -( ft + fc ) / ( 3. * fc ) / log(mQCompression / mQTension);
    const double mQ = mQTension * exp(R / AConst);

    //compute yieldHard and yieldSoft
    const double yieldHardOne = computeHardeningOne(tempKappa);

    const double Bl = sig / fc;

    const double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.);

    const double dgdsig = 4. * ( 1. - yieldHardOne ) / fc * Al * Bl + yieldHardOne * yieldHardOne * mQ / fc;

    const double dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +
                          m *pow(yieldHardOne, 2.) / ( sqrt(6.) * fc );

    /* Debug: Introduce the factor 2G/K. I am not sure if this is correct.
     * There might be too many stress states in the vertex case. */
    return dgdrho / dgdsig * 3. * ( 1. - 2. * nu ) / ( 1. + nu );
}


void
ConcreteDPM :: computeDDGDInvDKappa(FloatArray &answer,
                                    const double sig,
                                    const double rho,
                                    const double tempKappa)
{
    //Compute dilation parameter
    const double R = ( sig - ft / 3. ) / fc;
    const double mQTension = ( 3. * ft / fc + m / 2. );
    const double mQCompression =
        ( 1. + 2. * dilationConst ) / ( dilationConst - 1. ) * ( 3. + m / 2. );
    const double AConst = -( ft + fc ) / ( 3. * fc ) / log(mQCompression / mQTension);
    const double mQ = mQTension * exp(R / AConst);

    //compute yieldHard and yieldSoft
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double dYieldHardOneDKappa = computeHardeningOnePrime(tempKappa);

    const double Bl = sig / fc + rho / ( fc * sqrt(6.) );

    const double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fc;

    const double dAlDYieldHard = -pow(Bl, 2.);

    const double dDGDSigDKappa =
        ( -4. * Al * Bl / fc + 4. * ( 1 - yieldHardOne ) / fc * dAlDYieldHard * Bl +
          2. * yieldHardOne * mQ / fc ) * dYieldHardOneDKappa;

    const double dDGDRhoDKappa =
        ( dAlDYieldHard / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) -
         4. * Al / ( sqrt(6.) * fc ) * Bl + 2. * m * yieldHardOne / ( sqrt(6.) * fc ) ) * dYieldHardOneDKappa;

    answer(0) = dDGDSigDKappa;
    answer(1) = dDGDRhoDKappa;
}

void
ConcreteDPM :: computeDDGDDInv(FloatMatrix &answer,
                               const double sig,
                               const double rho,
                               const double tempKappa)
{
    //Compute dilation parameter
    const double R = ( sig - ft / 3. ) / fc;
    const double mQTension = ( 3. * ft / fc + m / 2. );
    const double mQCompression =
        ( 1. + 2. * dilationConst ) / ( dilationConst - 1. ) * ( 3. + m / 2. );
    const double AConst = -( ft + fc ) / ( 3. * fc ) / log(mQCompression / mQTension);
    //  const double mQ = mQTension*exp(R/AConst);
    const double dMQDSig = mQTension / ( AConst * fc ) * exp(R / AConst);

    //compute yieldHardOne and yieldSoft
    const double yieldHardOne = computeHardeningOne(tempKappa);

    //compute help parameter Al and Bl and the corresponding derivatives
    const double Bl = sig / fc + rho / ( fc * sqrt(6.) );

    const double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) +
                      sqrt(3. / 2.) * rho / fc;

    const double dAlDSig = 2. * ( 1. - yieldHardOne ) * Bl / fc;
    const double dBlDSig = 1. / fc;

    const double dAlDRho = 2. * ( 1. - yieldHardOne ) * Bl / ( fc * sqrt(6.) ) + sqrt(3. / 2.) / fc;
    const double dBlDRho = 1. / ( fc * sqrt(6.) );

    //compute second derivatives of g
    const double ddgddSig = 4. * ( 1. - yieldHardOne ) / fc * ( dAlDSig * Bl + Al * dBlDSig ) +
                            yieldHardOne * yieldHardOne * dMQDSig / fc;

    const double ddgddRho = dAlDRho / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +
                            Al * dBlDRho * 4. * ( 1. - yieldHardOne ) / ( sqrt(6.) * fc );

    const double ddgdSigdRho = 4. * ( 1. - yieldHardOne ) / fc * ( dAlDRho * Bl + Al * dBlDRho );

    const double ddgdRhodSig = dAlDSig / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. )
                               + Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * dBlDSig );

    answer(0, 0) = ddgddSig;
    answer(0, 1) = ddgdSigdRho;
    answer(1, 0) = ddgdRhodSig;
    answer(1, 1) = ddgddRho;
}

void
ConcreteDPM :: computeAMatrix(FloatMatrix &answer,
                              const double sig,
                              const double rho,
                              const double tempKappa)
{
    FloatMatrix aMatrixInverse(3, 3);
    FloatArray dDKappaDDeltaLambdaDInv(2);
    FloatMatrix dDGDDInv(2, 2);
    FloatArray dDGDInvDKappa(2);
    double dDKappaDDeltaLambdaDKappa;
    answer.zero();
    computeDDGDDInv(dDGDDInv, sig, rho, tempKappa);
    computeDDKappaDDeltaLambdaDInv(dDKappaDDeltaLambdaDInv, sig, rho,
                                   tempKappa);
    dDKappaDDeltaLambdaDKappa =
        computeDDKappaDDeltaLambdaDKappa(sig, rho, tempKappa);
    computeDDGDInvDKappa(dDGDInvDKappa, sig, rho, tempKappa);


    aMatrixInverse(0, 0) = 1. / ( kM ) + deltaLambda *dDGDDInv(0, 0);
    aMatrixInverse(0, 1) = deltaLambda * dDGDDInv(0, 1);
    aMatrixInverse(0, 2) = deltaLambda * dDGDInvDKappa(0);

    aMatrixInverse(1, 0) = deltaLambda * dDGDDInv(1, 0);
    aMatrixInverse(1, 1) = 1. / ( 2. * gM ) + deltaLambda *dDGDDInv(1, 1);
    aMatrixInverse(1, 2) = deltaLambda * dDGDInvDKappa(1);

    aMatrixInverse(2, 0) = deltaLambda * dDKappaDDeltaLambdaDInv(0);
    aMatrixInverse(2, 1) = deltaLambda * dDKappaDDeltaLambdaDInv(1);
    aMatrixInverse(2, 2) = -1. + deltaLambda * dDKappaDDeltaLambdaDKappa;

    answer.beInverseOf(aMatrixInverse);
}



double
ConcreteDPM :: computeHardeningOne(const double kappa) const
{
    if ( kappa <= 0. ) {
        return yieldHardInitial;
    } else if ( kappa > 0. && kappa < 1. ) {
        return yieldHardInitial + ( 1. - yieldHardInitial ) *
               kappa * ( pow(kappa, 2.) - 3. * kappa + 3. );
    } else {
        return 1.;
    }
}


double
ConcreteDPM :: computeHardeningOnePrime(const double kappa) const
{
    if ( kappa < 0. ) {
        return 3. * ( 1. - yieldHardInitial );
    } else if ( kappa >= 0. && kappa < 1. ) {
        return ( 1. - yieldHardInitial ) * ( 3. * pow(kappa, 2.) - 6. * kappa + 3. );
    } else {
        return 0.;
    }
}

void
ConcreteDPM :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                             MatResponseForm form,
                                             MatResponseMode mode,
                                             GaussPoint *gp,
                                             TimeStep *atTime)
{
    if ( gp->giveMaterialMode() == _3dMat || gp->giveMaterialMode() ==  _PlaneStrain || gp->giveMaterialMode() == _3dRotContinuum  ) {
        double omega = 0.;
        ConcreteDPMStatus *status = giveStatus(gp);
        if ( mode == ElasticStiffness ) {
            this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
        } else if ( mode == SecantStiffness || mode == TangentStiffness ) {
            omega = status->giveTempDamage();
            if ( omega > 0.9999 ) {
                omega = 0.9999;
            }

            this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
            answer.times(1. - omega);
        }
    }
}

void
ConcreteDPM :: computeTrialCoordinates(const StressVector &stress)
{
    StressVector deviatoricStress(matMode);
    effectiveStress.computeDeviatoricVolumetricSplit(deviatoricStress,
                                                     sig);
    rho = deviatoricStress.computeSecondCoordinate();
    thetaTrial = deviatoricStress.computeThirdCoordinate();
}


void
ConcreteDPM :: assignStateFlag(GaussPoint *gp)
{
    ConcreteDPMStatus *status = giveStatus(gp);
    //Get kappaD from status to define state later on
    damage = status->giveDamage();
    tempDamage = status->giveTempDamage();
    kappaP = status->giveKappaP();
    tempKappaP = status->giveTempKappaP();

    if ( tempKappaP > kappaP ) {
        if ( tempDamage > damage ) {
            status->
            letTempStateFlagBe(ConcreteDPMStatus :: ConcreteDPM_PlasticDamage);
        } else {
            status->letTempStateFlagBe(ConcreteDPMStatus :: ConcreteDPM_Plastic);
        }
    } else {
        const int state_flag = status->giveStateFlag();
        if ( state_flag == ConcreteDPMStatus :: ConcreteDPM_Elastic ) {
            if ( tempDamage > damage ) {
                status->
                letTempStateFlagBe(ConcreteDPMStatus :: ConcreteDPM_Damage);
            } else {
                status->
                letTempStateFlagBe(ConcreteDPMStatus :: ConcreteDPM_Elastic);
            }
        } else {
            if ( tempDamage > damage ) {
                status->
                letTempStateFlagBe(ConcreteDPMStatus :: ConcreteDPM_Damage);
            } else {
                status->
                letTempStateFlagBe(ConcreteDPMStatus :: ConcreteDPM_Unloading);
            }
        }
    }
}

void
ConcreteDPM :: computeDRhoDStress(FloatArray &answer,
                                  const StressVector &stress) const
{
    int size = 6;
    //compute volumetric deviatoric split
    StressVector deviatoricStress(_3dMat);
    double volumetricStress;
    stress.computeDeviatoricVolumetricSplit(deviatoricStress, volumetricStress);
    double rho = deviatoricStress.computeSecondCoordinate();

    //compute the derivative of J2 with respect to the stress
    FloatArray dJ2DStress;
    dJ2DStress = deviatoricStress;
    for ( int i = 3; i < size; i++ ) {
        dJ2DStress(i) = deviatoricStress(i) * 2.0;
    }

    //compute the derivative of rho with respect to stress
    FloatArray dRhoDStress;
    dRhoDStress = dJ2DStress;
    dRhoDStress.times(1. / rho);

    answer = dRhoDStress;
}

void
ConcreteDPM :: computeDSigDStress(FloatArray &answer) const
{
    int size = 6;
    for ( int i = 0; i < 3; i++ ) {
        answer(i) = 1. / 3.;
    }

    for ( int i = 3; i < size; i++ ) {
        answer(i) = 0.;
    }
}


void
ConcreteDPM :: computeDDRhoDDStress(FloatMatrix &answer,
                                    const StressVector &stress) const

{
    int size = 6;

    //compute volumetric deviatoric split
    StressVector deviatoricStress(_3dMat);
    double volumetricStress;
    stress.computeDeviatoricVolumetricSplit(deviatoricStress, volumetricStress);
    double rho = deviatoricStress.computeSecondCoordinate();


    //compute first derivative of J2
    FloatArray dJ2dstress;
    dJ2dstress = deviatoricStress;
    for ( int i = 3; i < deviatoricStress.giveSize(); i++ ) {
        dJ2dstress(i) = deviatoricStress(i) * 2.;
    }

    //compute second derivative of J2
    FloatMatrix ddJ2ddstress(size, size);
    ddJ2ddstress.zero();
    for ( int i = 0; i < size; i++ ) {
        if ( i < 3 ) {
            ddJ2ddstress(i, i) = 2. / 3.;
        }

        if ( i > 2 ) {
            ddJ2ddstress(i, i) = 2.;
        }
    }

    ddJ2ddstress(0, 1) = -1. / 3.;
    ddJ2ddstress(0, 2) = -1. / 3.;
    ddJ2ddstress(1, 0) = -1. / 3.;
    ddJ2ddstress(1, 2) = -1. / 3.;
    ddJ2ddstress(2, 0) = -1. / 3.;
    ddJ2ddstress(2, 1) = -1. / 3.;

    //compute square of the first derivative of J2
    FloatMatrix dJ2DJ2(size, size);
    for ( int v = 0; v < size; v++ ) {
        for ( int w = 0; w < size; ++w ) {
            dJ2DJ2(v, w) = dJ2dstress(v) * dJ2dstress(w);
        }
    }

    //compute the second derivative of rho
    FloatMatrix ddRhoddStress;
    ddRhoddStress = ddJ2ddstress;
    ddRhoddStress.times(1. / rho);
    FloatMatrix help1;
    help1 = dJ2DJ2;
    help1.times( -1. / ( rho * rho * rho ) );
    ddRhoddStress.add(help1);
    answer = ddRhoddStress;
}

void
ConcreteDPM :: computeDCosThetaDStress(FloatArray &answer,
                                       const StressVector &stress) const
{
    int size = stress.giveSize();

    //compute volumetric deviatoric split
    StressVector deviatoricStress(_3dMat);
    double volumetricStress;
    stress.computeDeviatoricVolumetricSplit(deviatoricStress, volumetricStress);

    //compute the coordinates
    double rho = deviatoricStress.computeSecondCoordinate();

    //compute principal stresses and directions
    FloatArray principalDeviatoricStress;
    FloatMatrix principalDir;
    deviatoricStress.computePrincipalValDir(principalDeviatoricStress, principalDir);

    //compute the derivative of s1 with respect to the cartesian stress
    FloatArray ds1DStress(size);
    ds1DStress(0) = principalDir(0, 0) * principalDir(0, 0) - 1. / 3.;
    ds1DStress(1) = principalDir(1, 0) * principalDir(1, 0) - 1. / 3.;
    ds1DStress(2) = principalDir(2, 0) * principalDir(2, 0) - 1. / 3.;
    ds1DStress(3) = 2. * principalDir(1, 0) * principalDir(2, 0);
    ds1DStress(4) = 2. * principalDir(2, 0) * principalDir(0, 0);
    ds1DStress(5) = 2. * principalDir(0, 0) * principalDir(1, 0);

    //compute dCosThetaDStress
    FloatArray dCosThetaDStress;
    dCosThetaDStress = ds1DStress;
    dCosThetaDStress.times( sqrt(3. / 2.) * rho / pow(rho, 2.) );
    FloatArray help(size);
    computeDRhoDStress(help, stress);
    help.times( -sqrt(3. / 2.) * principalDeviatoricStress(0) / pow(rho, 2.) );
    dCosThetaDStress.add(help);
    answer = dCosThetaDStress;
}

double
ConcreteDPM :: computeDRDCosTheta(const double theta, const double ecc) const
{
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

int
ConcreteDPM :: setIPValue(const FloatArray value, GaussPoint *gp, InternalStateType type)
{
    ConcreteDPMStatus *status = giveStatus(gp);
    if ( status->setIPValue(value, type) ) {
        return 1;
    } else {
        return StructuralMaterial :: setIPValue(value, gp, type);
    }
}

int
ConcreteDPM :: giveIPValue(FloatArray &answer,
                           GaussPoint *gp,
                           InternalStateType type,
                           TimeStep *atTime)
{
    const ConcreteDPMStatus *status = giveStatus(gp);
    //int state_flag = status->giveStateFlag();
    //double stateFlagValue = 0.;
    StrainVector plasticStrainVector(_3dMat);

    if ( type == IST_PlasticStrainTensor ) {
        answer.resize(6);
        status->giveFullPlasticStrainVector(plasticStrainVector);
        answer = plasticStrainVector;
        return 1;
    }

    if ( ( type == IST_DamageScalar ) || ( type == IST_DamageTensor ) || ( type == IST_PrincipalDamageTensor ) ) {
        answer.resize(1);
        answer.at(1) = status->giveDamage();
        return 1;
    }

    if ( ( type == IST_DamageTensorTemp ) || ( type == IST_PrincipalDamageTempTensor ) ) {
        answer.resize(1);
        answer.at(1) = status->giveTempDamage();
        return 1;
    }

    if ( type == IST_CumPlasticStrain ) {
        answer.resize(1);
        answer.at(1) = status->giveKappaP();
        return 1;
    }

    if ( type == IST_CumPlasticStrain_2 ) {
        answer.resize(1);
        answer.at(1) = status->giveKappaD();
        return 1;
    }

    if ( type == IST_VolumetricPlasticStrain ) {
        answer.resize(1);
        answer.at(1) = status->giveVolumetricPlasticStrain();
        return 1;
    }

    return StructuralMaterial :: giveIPValue(answer, gp, type, atTime);
}

int
ConcreteDPM :: giveIPValueSize(InternalStateType type,
                               GaussPoint *gp)
{
    if ( type == IST_PlasticStrainTensor ) {
      //return 6;       
	return this->giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() );
    } else if ( ( type == IST_CumPlasticStrain ) || ( type == IST_CumPlasticStrain_2 ) || ( type == IST_VolumetricPlasticStrain ) || ( type == IST_PrincipalDamageTensor ) || ( type == IST_PrincipalDamageTempTensor ) || ( type == IST_DamageScalar ) || ( type == IST_DamageTensor ) || ( type == IST_DamageTensorTemp ) ) {
        return 1;
    } else {
        return StructuralMaterial :: giveIPValueSize(type, gp);
    }
}

int
ConcreteDPM :: giveIntVarCompFullIndx(IntArray &answer,
                                      InternalStateType type,
                                      MaterialMode mmode)
{
    switch ( type ) {
    case IST_PlasticStrainTensor:
      this->giveStressStrainMask(answer, FullForm, mmode);
      /*
        answer.resize(6);
        answer.zero();
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 4;
        answer.at(5) = 5;
        answer.at(6) = 6;
      */
        return 1;

    case IST_DamageTensor:
    case IST_DamageTensorTemp:
        answer.resize(6);
        answer.zero();
        answer.at(1) = 1;
        return 1;

    case IST_DamageScalar:
    case IST_CumPlasticStrain:
    case IST_CumPlasticStrain_2:
    case IST_VolumetricPlasticStrain:
        answer.resize(1);
        answer.at(1) = 1;
        return 1;

    default:
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);

        break;
    }
}

InternalStateValueType
ConcreteDPM :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_PlasticStrainTensor ) {
        return ISVT_TENSOR_S3E;
    } else if ( ( type == IST_DamageTensor ) || ( type == IST_DamageTensorTemp ) ) {
        return ISVT_TENSOR_S3;
    } else if ( ( type == IST_DamageScalar ) || ( type == IST_CumPlasticStrain ) || ( type == IST_CumPlasticStrain_2 ) || ( type == IST_VolumetricPlasticStrain ) ) {
        return ISVT_SCALAR;
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }
}

MaterialStatus *
ConcreteDPM :: CreateStatus(GaussPoint *gp) const
{
    ConcreteDPMStatus *status =
        new  ConcreteDPMStatus(1, StructuralMaterial :: giveDomain(), gp);
    return status;
}
} // end namespace oofem
