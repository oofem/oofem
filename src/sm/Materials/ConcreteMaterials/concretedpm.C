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

#include "concretedpm.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "intarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(ConcreteDPM);
///@todo Eventually remove this old input string (replacing the name in input files is easy anyway).
//static bool __dummy_ConcreteDPM_alt __attribute__((unused)) = GiveClassFactory().registerMaterial("concreteidp", matCreator< ConcreteDPM > );
REGISTER_Material_Alt(ConcreteDPM, concreteidp);

ConcreteDPMStatus :: ConcreteDPMStatus(GaussPoint *gp) :
    StructuralMaterialStatus(gp)
{
    strainVector.resize(6);
    stressVector.resize(6);
    tempStressVector = stressVector;
    tempStrainVector = strainVector;
}


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
    tempEpsloc = epsloc;
}

void
ConcreteDPMStatus :: updateYourself(TimeStep *tStep)
{
#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    epsloc = tempEpsloc;
#endif
    // Call corresponding function of the parent class to update
    // variables defined there.
    StructuralMaterialStatus :: updateYourself(tStep);

    // update variables defined in ConcreteDPMStatus
    plasticStrain = tempPlasticStrain;
    kappaP = tempKappaP;
    kappaD = tempKappaD;
    damage = tempDamage;
    equivStrain = tempEquivStrain;
    state_flag = temp_state_flag;
    epsloc = tempEpsloc;
}

void
ConcreteDPMStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
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

#if 0
     // print plastic strain vector
     FloatArray plasticStrainVector();
     giveFullPlasticStrainVector(plasticStrainVector);
     fprintf(file,"plastic strains ");
     for ( auto &val : plasticStrainVector )
        fprintf(file," %.4e", val);
#endif

    // print hardening/softening parameters and damage

    //fprintf (file," kappaP %.4e", kappaP ) ;
    //fprintf (file,", kappaD %.4e", kappaD ) ;
    fprintf(file, ", kappa %.4e %.4e", kappaP, kappaD);
    fprintf(file, ", damage %.4e", damage);

    // end of record
    fprintf(file, "}\n");
}

void
ConcreteDPMStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = plasticStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(kappaP) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaD) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(equivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(state_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(deltaEquivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


void
ConcreteDPMStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus :: restoreContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = plasticStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(kappaP) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappaD) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(equivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(state_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(deltaEquivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}

int
ConcreteDPMStatus :: setIPValue(const FloatArray &value, InternalStateType type)
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


//   ******************************************************
//   *** CLASS CONCRETE DAMAGE-PLASTIC MATERIAL MODEL   ***
//   ******************************************************

//#define DPM_ITERATION_LIMIT 1.e-8

ConcreteDPM :: ConcreteDPM(int n, Domain *d) :
    StructuralMaterial(n, d),
    linearElasticMaterial(n,d)
{}


void
ConcreteDPM :: initializeFrom(InputRecord &ir)
{
    // call the corresponding service for the linear elastic material
    StructuralMaterial :: initializeFrom(ir);
    linearElasticMaterial.initializeFrom(ir);

    double value;
    // elastic parameters
    IR_GIVE_FIELD(ir, eM, _IFT_IsotropicLinearElasticMaterial_e);
    IR_GIVE_FIELD(ir, nu, _IFT_IsotropicLinearElasticMaterial_n);
    propertyDictionary.add('E', eM);
    propertyDictionary.add('n', nu);

    IR_GIVE_FIELD(ir, value, _IFT_IsotropicLinearElasticMaterial_talpha);
    propertyDictionary.add(tAlpha, value);

    gM = eM / ( 2. * ( 1. + nu ) );
    kM = eM / ( 3. * ( 1. - 2. * nu ) );

    // instanciate the variables of the plasticity model
    IR_GIVE_FIELD(ir, fc, _IFT_ConcreteDPM_fc);
    IR_GIVE_FIELD(ir, ft, _IFT_ConcreteDPM_ft);
    propertyDictionary.add(ft_strength, ft);
    propertyDictionary.add(fc_strength, fc);

    // damage parameters - only exponential softening
    // [in ef variable the wf (crack opening) is stored]
    if ( ir.hasField(_IFT_ConcreteDPM_wf) ) {
        IR_GIVE_FIELD(ir, ef, _IFT_ConcreteDPM_wf);
        // fracture energy
    } else {
        double Gf;
        IR_GIVE_FIELD(ir, Gf, _IFT_ConcreteDPM_gf);
        ef = Gf / ft; // convert fracture energy to crack opening
    }

    // default parameters
    ecc = 0.525;
    IR_GIVE_OPTIONAL_FIELD(ir, ecc, _IFT_ConcreteDPM_ecc);
    yieldHardInitial = 0.1;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldHardInitial, _IFT_ConcreteDPM_kinit);
    AHard = 8.e-2;
    IR_GIVE_OPTIONAL_FIELD(ir, AHard, _IFT_ConcreteDPM_ahard);
    BHard = 3.e-3;
    IR_GIVE_OPTIONAL_FIELD(ir, BHard, _IFT_ConcreteDPM_bhard);
    CHard = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, CHard, _IFT_ConcreteDPM_chard);
    DHard = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, DHard, _IFT_ConcreteDPM_dhard);
    ASoft = 15.;
    IR_GIVE_OPTIONAL_FIELD(ir, ASoft, _IFT_ConcreteDPM_asoft);
    helem = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, helem, _IFT_ConcreteDPM_helem);
#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    href = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, href, _IFT_ConcreteDPM_href);
#endif

    //Compute m
    m = 3. * ( pow(fc, 2.) - pow(ft, 2.) ) / ( fc * ft ) * ecc / ( ecc + 1. );
    //Compute default value of dilationConst
    dilationConst = -0.85;
    IR_GIVE_OPTIONAL_FIELD(ir, dilationConst, _IFT_ConcreteDPM_dilation);

    yieldTol = 1.e-10;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, _IFT_ConcreteDPM_yieldtol);
    newtonIter = 100;
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, _IFT_ConcreteDPM_newtoniter);
}

FloatArrayF<6>
ConcreteDPM :: giveRealStressVector_3d(const FloatArrayF<6> &totalStrain, GaussPoint *gp,
                                       TimeStep *tStep) const
{
    auto status = giveConcreteDPMStatus(gp);

    // Initialize temp variables for this gauss point
    status->initTempStatus();


    // subtract stress-independent part of strain
    // (due to temperature changes, shrinkage, etc.)
    auto thermalStrain = this->computeStressIndependentStrainVector_3d(gp, tStep, VM_Total);
    auto strain = totalStrain - thermalStrain;

    // perform plasticity return
    performPlasticityReturn(gp, strain);

    // compute damage
    //auto [tempDamage, tempKappaD] = computeDamage(strain, gp, tStep); // c++17
    auto tempDamKapD = computeDamage(strain, gp, tStep); // c++17
    auto tempDamage = tempDamKapD.first;
    auto tempKappaD = tempDamKapD.second;

    // compute elastic strains and effective stress
    const auto &tempPlasticStrain = status->giveTempPlasticStrain();
    auto elasticStrain = strain - tempPlasticStrain;
    auto effectiveStress = applyElasticStiffness(elasticStrain, eM, nu);

    // compute the nominal stress
    auto stress = effectiveStress * (1. - tempDamage);

    status->letTempKappaDBe(tempKappaD);
    status->letTempDamageBe(tempDamage);

    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(stress);

    assignStateFlag(gp);



#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    // check whether the second-order work is negative
    // (must be done !!!before!!! the update of strain and stress)
    if ( status->giveEpsLoc() < 0. ) {
        auto strainIncrement = totalStrain - status->giveStrainVector();
        auto stressIncrement = stress - FloatArrayF<6>(status->giveStressVector());
        int n = strainIncrement.giveSize();
        double work = strainIncrement.dotProduct(stressIncrement, n);
        //printf(" work : %g\n", work);
        if ( work < 0. ) {
            double E = this->give('E', gp);
            double ft = this->give(ft_strength, gp);
            double tmpEpsloc = tempKappaD + tempDamage * ft / E;   /// is this right? or was is supposed to use the old converged values?
            status->letTempEpslocBe(tmpEpsloc);
        }
    }

#endif
    return stress;
}


std::pair<double, double>
ConcreteDPM :: computeDamage(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = giveConcreteDPMStatus(gp);
    double equivStrain = computeEquivalentStrain(strain, gp, tStep);
    double f = equivStrain - status->giveKappaD();
    if ( f <= 0.0 ) {
        // damage does not grow
        double tempKappaD = status->giveKappaD();
        double tempDamage = status->giveDamage();
        return {tempDamage, tempKappaD};
    } else {
        // damage grows
        double tempKappaD = equivStrain;
#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
        if ( ( href <= 0. ) || ( status->giveEpsLoc() >= 0. ) ) {
            // evaluate and store the effective element size (if not known yet)
            this->initDamaged(tempKappaD, strain, gp);
        }

#else
        this->initDamaged(tempKappaD, strain, gp);
#endif
        double tempDamage = computeDamageParam(tempKappaD, gp);
        return {tempDamage, tempKappaD};
    }
}

double
ConcreteDPM :: computeEquivalentStrain(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    //The equivalent  strain is based on the volumetric plastic strain
    auto status = giveConcreteDPMStatus(gp);
    auto tempKappaP = status->giveTempKappaP(); /// FIXME
    auto kappaP = status->giveKappaP();  /// FIXME
    double equivStrain = status->giveEquivStrain();
    double deltaEquivStrain = 0.;
    double tempEquivStrain = 0.;
    if ( tempKappaP <= 1.0 || tempKappaP == kappaP ) {
        return equivStrain;
    } else if ( tempKappaP > 1.0 && tempKappaP != kappaP ) {
        const auto &plasticStrain = status->givePlasticStrain();
        const auto &tempPlasticStrain = status->giveTempPlasticStrain();

        double volumetricPlasticStrain = plasticStrain[0] + plasticStrain[1] + plasticStrain[2];
        double tempVolumetricPlasticStrain = tempPlasticStrain[0] + tempPlasticStrain[1] + tempPlasticStrain[2];
        if ( kappaP < 1.0 ) {
            //compute volumetric plastic strain at peak
            double peakVolumetricPlasticStrain = ( 1. - kappaP ) / ( tempKappaP - kappaP ) *
                                                 ( tempVolumetricPlasticStrain - volumetricPlasticStrain ) +
                                                 volumetricPlasticStrain;
            if ( peakVolumetricPlasticStrain < 0. ) {
                peakVolumetricPlasticStrain = 0.;
            }

            deltaEquivStrain = tempVolumetricPlasticStrain - peakVolumetricPlasticStrain;
            tempEquivStrain = deltaEquivStrain / computeDuctilityMeasureDamage(strain, gp);
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
    return tempEquivStrain;
}

double
ConcreteDPM :: computeInverseDamage(double dam, GaussPoint *gp) const
{
    auto status = giveConcreteDPMStatus(gp);
    double le = status->giveLe();
    if ( le == 0. ) {
        if ( helem > 0. ) {
            le = helem;
        } else {
            le = gp->giveElement()->computeMeanSize();
        }

        status->setLe(le);
    }

    double answer = -( this->ef / le ) * log(1. - dam) - dam * ft / eM;
    return answer;
}

#define DPM_DAMAGE_TOLERANCE 1.e-8

double
ConcreteDPM :: computeDamageParam(double kappa, GaussPoint *gp) const
{
    double omega = 0.;
    if ( kappa > 0. ) {
        // iteration to achieve mesh objectivity
        // this is only valid for tension
        ConcreteDPMStatus *status = giveConcreteDPMStatus(gp);
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
            OOFEM_ERROR("element too large");
        }

        do {
            nite++;
            aux = exp(-h * ( omega * ft / eM + kappa ) / ef);
            R = 1. - omega - aux;
            Lhs = -1. + aux1 * aux;
            omega -= R / Lhs;
            if ( nite > 40 ) {
                OOFEM_ERROR("algorithm not converging");
            }
        } while ( fabs(R) >= DPM_DAMAGE_TOLERANCE );

#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    } else {       // after localization, more complicated formula
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
        } while ( fabs(R) >= DPM_DAMAGE_TOLERANCE );
    }

#endif
        if ( ( omega > 1.0 ) || ( omega < 0.0 ) ) {
            OOFEM_ERROR("internal error, omega = %g", omega);
        }
    }

    return omega;
}


void
ConcreteDPM :: initDamaged(double kappaD,
                           const FloatArrayF<6> &strain,
                           GaussPoint *gp) const
{
    auto status = giveConcreteDPMStatus(gp);

    if ( kappaD <= 0. ) {
        return;
    }

    if ( helem > 0. ) {
        status->setLe(helem);
    } else if ( status->giveDamage() == 0. ) {
        //auto [principalStrains, principalDir] = computePrincipalValDir(from_voigt_strain(strain)); // c++17
        auto tmp = computePrincipalValDir(from_voigt_strain(strain));
        auto principalStrains = tmp.first;
        auto principalDir = tmp.second;
        // find index of max positive principal strain
        int indx = 1;
        for ( int i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) > principalStrains.at(indx) ) {
                indx = i;
            }
        }

        FloatArrayF<3> crackPlaneNormal;
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
        status->setLe(gp->giveElement()->computeMeanSize());
    }
}

double
ConcreteDPM :: computeDuctilityMeasureDamage(const FloatArrayF<6> &strain, GaussPoint *gp) const
{
    auto status = giveConcreteDPMStatus(gp);
    const auto &plasticStrain = status->givePlasticStrain();
    auto tempPlasticStrain = status->giveTempPlasticStrain() - plasticStrain;

    double volStrain = tempPlasticStrain[0] + tempPlasticStrain[1] + tempPlasticStrain[2];
    auto principalStrain = computePrincipalValues(from_voigt_strain(strain));
    
    double negativeVolStrain = 0.;
    for ( int i = 0; i < 3; i++ ) {
        if ( principalStrain[i] < 0. ) {
            negativeVolStrain += principalStrain[i];
        }
    }

    //compute ductility measure
    double Rs = -negativeVolStrain / volStrain;
    if ( Rs < 1.0 ) {
        return 1. + ASoft *pow(Rs, 2.);
    } else {
        return 1. - 3. * ASoft + 4. *ASoft *sqrt(Rs);
    }
}

void
ConcreteDPM :: performPlasticityReturn(GaussPoint *gp,
                                       FloatArrayF<6> &strain) const
{
    auto status = static_cast< ConcreteDPMStatus * >( this->giveStatus(gp) );

    status->initTempStatus();

    //get temp plastic strain and tempKappa
    auto tempPlasticStrain = status->giveTempPlasticStrain();
    auto tempKappaP = status->giveTempKappaP(); /// FIXME

    // compute elastic strains and trial stress
    auto elasticStrain = strain - tempPlasticStrain;
    auto effectiveStress = applyElasticStiffness(elasticStrain, eM, nu); /// FIXME

    //Compute trial coordinates
    //auto [sig, rho, thetaTrial] = computeTrialCoordinates(effectiveStress, gp); // c++17
    auto tmp = computeTrialCoordinates(effectiveStress, gp); // c++17
    double sig = std::get<0>(tmp);
    double rho = std::get<1>(tmp);
    double thetaTrial = std::get<2>(tmp);

    double yieldValue = computeYieldValue(sig, rho, thetaTrial, tempKappaP);
    // choose correct stress return and update state flag
    if ( yieldValue  > yieldTol ) {
        double apexStress = 0.;
        auto vertexType = checkForVertexCase(apexStress, sig, tempKappaP);

        //Make the appropriate return
        if ( vertexType == VT_Tension || vertexType == VT_Compression ) {
            performVertexReturn(effectiveStress, apexStress, gp);
            if ( vertexType == VT_Regular ) {
                //This was no real vertex case
                //get the original tempKappaP and stress
                tempKappaP = status->giveTempKappaP();
                effectiveStress = applyElasticStiffness(elasticStrain, eM, nu);
            }
        }

        if ( vertexType == VT_Regular ) {
            performRegularReturn(effectiveStress, gp);
        }
    }

    // update temp kappaP
    status->letTempKappaPBe(tempKappaP);
    // compute the plastic strains
    elasticStrain = applyElasticCompliance(effectiveStress, eM, nu);
    tempPlasticStrain = strain - elasticStrain;
    status->letTempPlasticStrainBe(tempPlasticStrain);
}


ConcreteDPM :: Concrete_VertexType
ConcreteDPM :: checkForVertexCase(double &answer,
                                  const double sig,
                                  const double tempKappa) const
{
    //Compute sigZero and compare with actual sig
    double apexCompression = 0.;

    //compressive apex
    if ( ( tempKappa < 1. ) && ( sig < 0. ) ) {
        const double yieldHardOne = computeHardeningOne(tempKappa);
        double sigZero = -15. * fc;
        double FZero = 1.;
        int l = 0;
        while ( fabs(FZero) > yieldTol && l <= newtonIter ) {
            l++;
            FZero = pow( ( 1. - yieldHardOne ), 2. ) * pow( ( sigZero / fc ), 4. ) +
            pow(yieldHardOne, 2.) * m * ( sigZero / fc ) - pow(yieldHardOne, 2.);

            double dFZeroDSigZero = pow( ( 1. - yieldHardOne ), 2. ) * 4. * pow( ( sigZero / fc ), 3. ) / fc +
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
        answer = 0.;
        return VT_Tension;
    } else if ( sig < apexCompression ) {
        answer = apexCompression;
        return VT_Compression;
    } else {
        return VT_Regular;
    }
}

void
ConcreteDPM :: performRegularReturn(FloatArrayF<6> &effectiveStress, GaussPoint *gp) const
{
    auto status = static_cast< ConcreteDPMStatus * >( this->giveStatus(gp) );

    //compute the principal directions of the stress
    //auto [help, stressPrincipalDir] = computePrincipalValDir(from_voigt_stress(effectiveStress)); // c++17
    auto tmp = computePrincipalValDir(from_voigt_stress(effectiveStress));
    auto stressPrincipalDir = tmp.second;

    //compute invariants from stress state
    //auto [deviatoricStress, sig] = this->computeDeviatoricVolumetricSplit(effectiveStress); // c++17
    auto tmp2 = this->computeDeviatoricVolumetricSplit(effectiveStress);
    auto deviatoricStress = tmp2.first;
    double sig = tmp2.second; /// FIXME: verify
    double rho = computeSecondCoordinate(deviatoricStress); /// FIXME: verify
    double deltaLambda = 0.; /// FIXME: verify

    const double volumetricPlasticStrain = 3. * status->giveVolumetricPlasticStrain();

    const double deviatoricPlasticStrainNorm = status->giveDeviatoricPlasticStrainNorm();

    double tempVolumetricPlasticStrain = volumetricPlasticStrain;
    double tempDeviatoricPlasticStrainNorm = deviatoricPlasticStrainNorm;

    const double kappaP = status->giveKappaP();
    auto tempKappaP = kappaP; /// FIXME

    auto yieldValue = computeYieldValue(sig, rho, thetaTrial, tempKappaP);
    double residualNorm = fabs(yieldValue);

    int negativeRhoFlag = 0;
    
    int iterationCount = 0;
    while ( residualNorm > yieldTol ) {
        if ( ++iterationCount == newtonIter ) {
            OOFEM_ERROR("Closest point projection did not converge.");
        }

        //compute the stress, yield value and resdiuals
        auto dGDInv = computeDGDInv(sig, rho, tempKappaP);
        auto dFDInv = computeDFDInv(sig, rho, tempKappaP);
        auto dFDKappa = computeDFDKappa(sig, rho, tempKappaP);
        auto dKappaDDeltaLambda = computeDKappaDDeltaLambda(sig, rho, tempKappaP);
        yieldValue = computeYieldValue(sig, rho, thetaTrial, tempKappaP);
        //The residual volumetric plastic strain is defined as eps1+eps2+eps3
        FloatArrayF<3> residual;
        residual[0] = volumetricPlasticStrain - tempVolumetricPlasticStrain +
                      deltaLambda * dGDInv[0];
        residual[1] = deviatoricPlasticStrainNorm -
                      tempDeviatoricPlasticStrainNorm + deltaLambda * dGDInv[1];
        residual[2] = kappaP - tempKappaP + deltaLambda * dKappaDDeltaLambda;

        // weighted norm
        residualNorm = norm(residual);
        //    printf("\n residualNorm= %e\n", residualNorm);
        if ( residualNorm > yieldTol ) {
            auto aMatrix = computeAMatrix(sig, rho, tempKappaP, deltaLambda);

            // assemble the derivatives of the yield surface
            FloatArrayF<3> derivativesOfYieldSurface;
            for ( int i = 0; i < 2; i++ ) {
                derivativesOfYieldSurface[i] = dFDInv[i];
            }

            derivativesOfYieldSurface[2] = dFDKappa;
            //assemble flow rules
            FloatArrayF<3> flowRules;
            for ( int i = 0; i < 2; i++ ) {
                flowRules[i] = dGDInv[i];
            }

            flowRules[2] = dKappaDDeltaLambda;

            //compute the deltaLambdaIncrement
            double deltaLambdaIncrement = yieldValue;
            auto helpVectorA = dot(aMatrix, residual);
            for ( int i = 0; i < 3; i++ ) {
                deltaLambdaIncrement -= derivativesOfYieldSurface[i] * helpVectorA[i];
            }

            auto helpVectorB = dot(aMatrix, flowRules);
            double denominator = 0.;
            for ( int i = 0; i < 3; i++ ) {
                denominator += derivativesOfYieldSurface[i] * helpVectorB[i];
            }

            deltaLambdaIncrement /= denominator;

            //compute increment of answer
            auto answerIncrement = -dot(aMatrix, (residual + deltaLambdaIncrement * flowRules));
            auto rhoTest = rho + answerIncrement[1];

            //Special case, if rho changes sign
            double tempKappaPTest = 0.;
            if ( rhoTest < 0. && negativeRhoFlag == 0 ) {
                //Determine deltaLambdaIncrement, so that rho is equal to zero
                answerIncrement[1] = -rho;

                double deltaLambdaIncrementNew =
                    ( -aMatrix(1, 0) * residual[0] - aMatrix(1, 1) * residual[1] -
                     aMatrix(1, 2) * residual[2] - answerIncrement[1] ) /
                    ( flowRules[0] * aMatrix(1, 0) + flowRules[1] * aMatrix(1, 1) +
                     flowRules[2] * aMatrix(1, 2) );

                // Special case, if deltaLambdaIncrement is equal to zero.
                if ( fabs(deltaLambdaIncrementNew) < yieldTol * 1.e3 ) {
                    negativeRhoFlag = 1;
                    deltaLambdaIncrementNew = deltaLambdaIncrement;
                }

                answerIncrement = - dot(aMatrix, (residual + deltaLambdaIncrementNew * flowRules));

                sig += answerIncrement[0];
                rho += answerIncrement[1];

                tempKappaPTest = tempKappaP;
                tempKappaP += answerIncrement[2];
                deltaLambda += deltaLambdaIncrementNew;
            } else {
                tempKappaPTest = tempKappaP;
                tempKappaP += answerIncrement[2];
                sig += answerIncrement[0];
                rho += answerIncrement[1];

                deltaLambda += deltaLambdaIncrement;
            }

            //Special case, if deltaKappaP is negative
            if ( ( tempKappaP - status->giveKappaP() ) < 0. ) {
                tempKappaP = tempKappaPTest;
            }

            tempVolumetricPlasticStrain -= answerIncrement[0] / ( kM );
            tempDeviatoricPlasticStrainNorm -= answerIncrement[1] / ( 2. * gM );
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
    FloatArrayF<6> stressPrincipal;
    stressPrincipal[0] = sig + sqrt(2. / 3.) * rho * cos(thetaTrial);
    stressPrincipal[1] = sig + sqrt(2. / 3.) * rho * cos(thetaTrial - 2. * M_PI / 3.);
    stressPrincipal[2] = sig + sqrt(2. / 3.) * rho * cos(thetaTrial + 2. * M_PI / 3.);

    effectiveStress = transformStressVectorTo(stressPrincipalDir, stressPrincipal, 1);
}

void
ConcreteDPM :: performVertexReturn(FloatArrayF<6> &effectiveStress,
                                   double apexStress,
                                   GaussPoint *gp) const
{
    auto status = static_cast< ConcreteDPMStatus * >( this->giveStatus(gp) );
    double sig2 = 0.;

    //auto [deviatoricStressTrial, sigTrial] = this->computeDeviatoricVolumetricSplit(effectiveStress); // c++17
    auto tmp = this->computeDeviatoricVolumetricSplit(effectiveStress);
    auto deviatoricStressTrial = tmp.first;
    auto sigTrial = tmp.second;
    
    const double rhoTrial = computeSecondCoordinate(deviatoricStressTrial);

    double tempKappaP = status->giveTempKappaP();  /// FIXME verify
    const double kappaInitial = tempKappaP;

    if ( vertexType == VT_Tension ) {
        sig2 = -0.1 * ft;
    } else if ( vertexType == VT_Compression ) {
        sig2 = apexStress;
    }

    tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sigTrial);

    double yieldValue = computeYieldValue(sigTrial, 0., 0., tempKappaP);

    tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sig2);

    double yieldValueMid = computeYieldValue(sig2, 0., 0., tempKappaP);

    if ( yieldValue * yieldValueMid >= 0. ) {
        vertexType = VT_Regular;
        return;
    }

    double sigAnswer;
    double ratioPotential;

    double dSig;
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
            for ( int i = 0; i < 3; i++ ) {
                effectiveStress[i] = sigAnswer;
            }

            for ( int i = 3; i < effectiveStress.giveSize(); i++ ) {
                effectiveStress[i] = 0.;
            }

            ratioPotential = computeRatioPotential(sigAnswer, tempKappaP);

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
                                const double sig) const
{
    //This function is called, if stress state is in vertex case
    double equivalentDeltaPlasticStrain = sqrt( 1. / 9. * pow( ( sigTrial - sig ) / ( kM ), 2. ) + pow(rhoTrial / ( 2. * gM ), 2.) );

    double rho = 0.; /// FIXME: verify
    double thetaVertex = 0.;
    double ductilityMeasure = computeDuctilityMeasure(sig, rho, thetaVertex);

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
                               const double tempKappa) const
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

FloatArrayF<2>
ConcreteDPM :: computeDFDInv(const double sig,
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

    return {dfdsig, dfdrho};
}

double
ConcreteDPM :: computeDKappaDDeltaLambda(const double sig,
                                         const double rho,
                                         const double tempKappa) const
{
    //Variables
    auto dGDInv = computeDGDInv(sig, rho, tempKappa);

    double equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv[0], 2.) +
                               pow(dGDInv[1], 2.) );

    double ductilityMeasure = computeDuctilityMeasure(sig, rho, this->thetaTrial);
    double dKappaDDeltaLambda = equivalentDGDStress / ductilityMeasure;
    return dKappaDDeltaLambda;
}


FloatArrayF<2>
ConcreteDPM :: computeDDKappaDDeltaLambdaDInv(const double sig,
                                              const double rho,
                                              const double tempKappa) const
{
    //Compute first and second derivative of plastic potential
    auto dGDInv = computeDGDInv(sig, rho, tempKappa);
    auto dDGDDInv = computeDDGDDInv(sig, rho, tempKappa);

    //Compute equivalentDGDStress
    double equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv[0], 2.) +
                               pow(dGDInv[1], 2.) );

    //computeDuctilityMeasure
    double ductilityMeasure = computeDuctilityMeasure(sig, rho, this->thetaTrial);

    //Compute dEquivalentDGDStressDInv
    FloatArrayF<2> dEquivalentDGDStressDInv;
    dEquivalentDGDStressDInv[0] =
        ( 2. / 3. * dGDInv[0] * dDGDDInv(0, 0) + 2. * dGDInv[1] * dDGDDInv(1, 0) ) / ( 2. * equivalentDGDStress );
    dEquivalentDGDStressDInv[1] =
        ( 2. / 3. * dGDInv[0] * dDGDDInv(0, 1) + 2. * dGDInv[1] * dDGDDInv(1, 1) ) / ( 2. * equivalentDGDStress );

    // compute the derivative of
    auto dDuctilityMeasureDInv = computeDDuctilityMeasureDInv(sig, rho, tempKappa);

    FloatArrayF<2> answer;
    answer[0] = ( dEquivalentDGDStressDInv[0] * ductilityMeasure - equivalentDGDStress * dDuctilityMeasureDInv[0] ) / pow(ductilityMeasure, 2.);
    answer[1] = ( dEquivalentDGDStressDInv[1] * ductilityMeasure - equivalentDGDStress * dDuctilityMeasureDInv[1] ) / pow(ductilityMeasure, 2.);
    return answer;
}

double
ConcreteDPM :: computeDDKappaDDeltaLambdaDKappa(const double sig,
                                                const double rho,
                                                const double tempKappa) const
{
    //Compute first and second derivative of plastic potential
    auto dGDInv = computeDGDInv(sig, rho, tempKappa);
    auto dDGDInvDKappa = computeDDGDInvDKappa(sig, rho, tempKappa);

    double equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv[0], 2.) + pow(dGDInv[1], 2.) );

    //computeDuctilityMeasure
    double ductilityMeasure = computeDuctilityMeasure(sig, rho, this->thetaTrial);

    //Compute dEquivalentDGDStressDKappa
    double dEquivalentDGDStressDKappa =
        ( 2. / 3. * dGDInv[0] * dDGDInvDKappa[0] + 2. * dGDInv[1] * dDGDInvDKappa[1] ) / ( 2. * equivalentDGDStress );

    // compute the derivative of
    double dDuctilityMeasureDKappa = 0.;

    double dDKappaDDeltaLambdaDKappa =
        ( dEquivalentDGDStressDKappa * ductilityMeasure -
          equivalentDGDStress * dDuctilityMeasureDKappa ) / pow(ductilityMeasure, 2.);

    return dDKappaDDeltaLambdaDKappa;
}

double
ConcreteDPM :: computeDuctilityMeasure(const double sig,
                                       const double rho,
                                       const double theta) const
{
    double thetaConst = pow(2. * cos(theta), 2.);
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

FloatArrayF<2>
ConcreteDPM :: computeDDuctilityMeasureDInv(const double sig,
                                            const double rho,
                                            const double tempKappa) const
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
        return {dDuctilityMeasureDX * dXDSig, 0.};
    } else {
        double dXDSig = -1. / fc;
        double dDuctilityMeasureDX = -( BHard - AHard ) / ( CHard ) / thetaConst *exp( -x / ( CHard ) );
        return {dDuctilityMeasureDX * dXDSig, 0.};
    }
}


FloatArrayF<2>
ConcreteDPM :: computeDGDInv(const double sig,
                             const double rho,
                             const double tempKappa) const
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

    return {dgdsig, dgdrho};
}


double
ConcreteDPM :: computeRatioPotential(const double sig,
                                     const double tempKappa) const
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


FloatArrayF<2>
ConcreteDPM :: computeDDGDInvDKappa(const double sig,
                                    const double rho,
                                    const double tempKappa) const
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

    return {dDGDSigDKappa, dDGDRhoDKappa};
}

FloatMatrixF<2,2>
ConcreteDPM :: computeDDGDDInv(const double sig,
                               const double rho,
                               const double tempKappa) const
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

    FloatMatrixF<2,2> answer;
    answer(0, 0) = ddgddSig;
    answer(0, 1) = ddgdSigdRho;
    answer(1, 0) = ddgdRhodSig;
    answer(1, 1) = ddgddRho;
    return answer;
}

FloatMatrixF<3,3>
ConcreteDPM :: computeAMatrix(const double sig,
                              const double rho,
                              const double tempKappa,
                              const double deltaLambda) const
{
    auto dDGDDInv = computeDDGDDInv(sig, rho, tempKappa);
    auto dDKappaDDeltaLambdaDInv = computeDDKappaDDeltaLambdaDInv(sig, rho, tempKappa);
    auto dDKappaDDeltaLambdaDKappa = computeDDKappaDDeltaLambdaDKappa(sig, rho, tempKappa);
    auto dDGDInvDKappa = computeDDGDInvDKappa(sig, rho, tempKappa);

    FloatMatrixF<3,3> aMatrixInverse;
    aMatrixInverse(0, 0) = 1. / ( kM ) + deltaLambda * dDGDDInv(0, 0);
    aMatrixInverse(0, 1) = deltaLambda * dDGDDInv(0, 1);
    aMatrixInverse(0, 2) = deltaLambda * dDGDInvDKappa[0];

    aMatrixInverse(1, 0) = deltaLambda * dDGDDInv(1, 0);
    aMatrixInverse(1, 1) = 1. / ( 2. * gM ) + deltaLambda * dDGDDInv(1, 1);
    aMatrixInverse(1, 2) = deltaLambda * dDGDInvDKappa[1];

    aMatrixInverse(2, 0) = deltaLambda * dDKappaDDeltaLambdaDInv[0];
    aMatrixInverse(2, 1) = deltaLambda * dDKappaDDeltaLambdaDInv[1];
    aMatrixInverse(2, 2) = -1. + deltaLambda * dDKappaDDeltaLambdaDKappa;

    return inv(aMatrixInverse);
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

FloatMatrixF<6,6>
ConcreteDPM :: give3dMaterialStiffnessMatrix(MatResponseMode mode,
                                             GaussPoint *gp,
                                             TimeStep *tStep) const
{
    auto status = giveConcreteDPMStatus(gp);
    auto d = this->linearElasticMaterial.give3dMaterialStiffnessMatrix(mode, gp, tStep);
    if ( mode == SecantStiffness || mode == TangentStiffness ) {
        double omega = status->giveTempDamage();
        if ( omega > 0.9999 ) {
            omega = 0.9999;
        }

        return d * (1. - omega);
    } else {
        return d;
    }
}
 

std::tuple<double, double, double>
ConcreteDPM :: computeTrialCoordinates(const FloatArrayF<6> &stress, GaussPoint *gp) const
{
    //auto [deviatoricStress, sig] = this->computeDeviatoricVolumetricSplit(stress); // c++17
    auto tmp = this->computeDeviatoricVolumetricSplit(stress);
    auto deviatoricStress = tmp.first;
    double sig = tmp.second;
    double rho = computeSecondCoordinate(deviatoricStress);
    double thetaTrial = computeThirdCoordinate(deviatoricStress);
    return {sig, rho, thetaTrial};
}


void
ConcreteDPM :: assignStateFlag(GaussPoint *gp) const
{
    auto status = giveConcreteDPMStatus(gp);
    //Get kappaD from status to define state later on
    const auto &damage = status->giveDamage();
    const auto &tempDamage = status->giveTempDamage();
    const auto &kappaP = status->giveKappaP();
    const auto &tempKappaP = status->giveTempKappaP();

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

FloatArrayF<6>
ConcreteDPM :: computeDRhoDStress(const FloatArrayF<6> &stress) const
{
    //compute volumetric deviatoric split
    auto deviatoricStress = this->computeDeviator(stress);
    double rho = computeSecondCoordinate(deviatoricStress);

    //compute the derivative of J2 with respect to the stress
    auto dJ2DStress = deviatoricStress;
    for ( int i = 3; i < 6; i++ ) {
        dJ2DStress[i] = deviatoricStress[i] * 2.0;
    }

    //compute the derivative of rho with respect to stress
    auto dRhoDStress = dJ2DStress * (1. / rho);
    return dRhoDStress;
}

void
ConcreteDPM :: computeDSigDStress(FloatArrayF<6> &answer) const
{
    int size = 6;
    for ( int i = 0; i < 3; i++ ) {
        answer[i] = 1. / 3.;
    }

    for ( int i = 3; i < size; i++ ) {
        answer[i] = 0.;
    }
}


FloatMatrixF<6,6>
ConcreteDPM :: computeDDRhoDDStress(const FloatArrayF<6> &stress) const
{
    //compute volumetric deviatoric split
    auto deviatoricStress = this->computeDeviator(stress);
    double rho = computeSecondCoordinate(deviatoricStress);

    //compute first derivative of J2
    auto dJ2dstress = deviatoricStress;
    for ( int i = 3; i < 6; i++ ) {
        dJ2dstress[i] = deviatoricStress[i] * 2.;
    }

    //compute second derivative of J2
    FloatMatrixF<6,6> ddJ2ddstress;
    for ( int i = 0; i < 6; i++ ) {
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
    auto dJ2DJ2 = dyad(dJ2dstress, dJ2dstress);

    //compute the second derivative of rho
    auto ddRhoddStress = ddJ2ddstress * (1. / rho) + dJ2DJ2 * ( -1. / ( rho * rho * rho ) );
    return ddRhoddStress;
}

FloatArrayF<6>
ConcreteDPM :: computeDCosThetaDStress(const FloatArrayF<6> &stress) const
{
    //compute volumetric deviatoric split
    auto deviatoricStress = computeDeviator(stress);

    //compute the coordinates
    double rho = computeSecondCoordinate(deviatoricStress);

    //compute principal stresses and directions
    //auto [principalDeviatoricStress, principalDir] = computePrincipalValDir(from_voigt_stress(deviatoricStress)); // c++17
    auto tmp = computePrincipalValDir(from_voigt_stress(deviatoricStress));
    auto principalDeviatoricStress = tmp.first;
    auto principalDir = tmp.second;

    //compute the derivative of s1 with respect to the cartesian stress
    FloatArrayF<6> ds1DStress;
    ds1DStress[0] = principalDir(0, 0) * principalDir(0, 0) - 1. / 3.;
    ds1DStress[1] = principalDir(1, 0) * principalDir(1, 0) - 1. / 3.;
    ds1DStress[2] = principalDir(2, 0) * principalDir(2, 0) - 1. / 3.;
    ds1DStress[3] = 2. * principalDir(1, 0) * principalDir(2, 0);
    ds1DStress[4] = 2. * principalDir(2, 0) * principalDir(0, 0);
    ds1DStress[5] = 2. * principalDir(0, 0) * principalDir(1, 0);

    //compute dCosThetaDStress
    auto dCosThetaDStress = ds1DStress * ( sqrt(3. / 2.) * rho / pow(rho, 2.) ) + 
        computeDRhoDStress(stress) * ( -sqrt(3. / 2.) * principalDeviatoricStress[0] / pow(rho, 2.) );
    return dCosThetaDStress;
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


void
ConcreteDPM :: restoreConsistency(GaussPoint *gp)
{
    auto status = giveConcreteDPMStatus(gp);

    // compute kappaD from damage
    double kappaD = this->computeInverseDamage(status->giveDamage(), gp);
    status->letKappaDBe(kappaD);
    status->letEquivStrainBe(kappaD);

    const auto &damage = status->giveDamage();

    // compute plastic strain
    // such that the given stress is obtained at zero total strain and given damage
    if ( damage < 1. ) {
        FloatArrayF<6> effectiveStress = status->giveStressVector();
        effectiveStress *= -1. / ( 1. - damage );
        auto D = this->give3dMaterialStiffnessMatrix(ElasticStiffness, gp, nullptr);
        auto plasticStrain = solve(D, effectiveStress);
        status->letPlasticStrainBe(plasticStrain);
    }
}


int
ConcreteDPM :: setIPValue(const FloatArray &value, GaussPoint *gp, InternalStateType type)
{
    auto status = giveConcreteDPMStatus(gp);
    if ( status->setIPValue(value, type) ) {
        return 1;
    } else {
        return StructuralMaterial :: setIPValue(value, gp, type);
    }
}

int
ConcreteDPM :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    const auto status = giveConcreteDPMStatus(gp);
    //int state_flag = status->giveStateFlag();
    //double stateFlagValue = 0.;

    if ( type == IST_PlasticStrainTensor ) {
        answer = status->givePlasticStrain();
        return 1;
    } else if ( type == IST_DamageScalar ) {
        answer.resize(1);
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
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveDamage();
        return 1;
    } else if ( type == IST_DamageTensorTemp ) {
        answer.resize(1);
        answer.at(1) = status->giveTempDamage();
        return 1;
    } else if ( type == IST_CumPlasticStrain ) {
        answer.resize(1);
        answer.at(1) = status->giveKappaP();
        return 1;
    } else if ( type == IST_CumPlasticStrain_2 ) {
        answer.resize(1);
        answer.at(1) = status->giveKappaD();
        return 1;
    } else if ( type == IST_VolumetricPlasticStrain ) {
        answer.resize(1);
        answer.at(1) = status->giveVolumetricPlasticStrain();
        return 1;
    }

    return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
}

MaterialStatus *
ConcreteDPM :: CreateStatus(GaussPoint *gp) const
{
    return new ConcreteDPMStatus(gp);
}
} // end namespace oofem
