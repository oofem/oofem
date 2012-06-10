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

#include "rcsd.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "isolinearelasticmaterial.h"
#include "datastream.h"
#include "contextioerr.h"

#ifndef __MAKEDEPEND
 #include <cstring>
#endif

namespace oofem {
RCSDMaterial :: RCSDMaterial(int n, Domain *d) : RCM2Material(n, d)
//
// constructor
//
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
}


RCSDMaterial :: ~RCSDMaterial()
//
// destructor
//
{
    delete linearElasticMaterial;
}


void
RCSDMaterial :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                     const FloatArray &totalStrain,
                                     TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    int i;
    FloatMatrix Ds0;
    double equivStrain;
    FloatArray princStress, crackStrain, reducedAnswer;
    FloatArray reducedStrainVector, strainVector, principalStrain;
    FloatArray reducedSpaceStressVector;
    FloatMatrix tempCrackDirs;
    RCSDMaterialStatus *status = ( RCSDMaterialStatus * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedStrainVector, gp, totalStrain,
                                                atTime, VM_Total);
    //

    crossSection->giveFullCharacteristicVector(strainVector, gp, reducedStrainVector);

    status->giveTempCrackDirs(tempCrackDirs);
    this->computePrincipalValDir(principalStrain, tempCrackDirs,
                                 strainVector,
                                 principal_strain);


    if ( status->giveTempMode() == RCSDMaterialStatus :: rcMode ) {
        // rotating crack mode

        this->giveRealPrincipalStressVector3d(princStress, gp, principalStrain, tempCrackDirs, atTime);
        princStress.resize(6);
        status->giveTempCrackDirs(tempCrackDirs);
        this->transformStressVectorTo(answer, tempCrackDirs, princStress, 1);

        crossSection->giveReducedCharacteristicVector(reducedSpaceStressVector, gp, answer);
        status->letTempStressVectorBe(reducedSpaceStressVector);

        status->giveCrackStrainVector(crackStrain);
        this->updateCrackStatus(gp, crackStrain);

        if ( form == ReducedForm ) {
            crossSection->giveReducedCharacteristicVector(reducedAnswer, gp, answer);
            answer = reducedAnswer;
        }

        // test if transition to scalar damage mode take place
        double minSofteningPrincStress = this->Ft, dCoeff, CurrFt, E, ep, ef, damage;
        int ipos = 0;
        for ( i = 1; i <= 3; i++ ) {
            if ( status->giveTempCrackStatus(i) == pscm_SOFTENING ) {
                if ( princStress.at(i) < minSofteningPrincStress ) {
                    minSofteningPrincStress = princStress.at(i);
                    ipos = i;
                }
            }
        }

        CurrFt = this->computeStrength( gp, status->giveCharLength(ipos) );

        if ( minSofteningPrincStress <= this->SDTransitionCoeff * CurrFt ) {
            // sd transition takes place

            E = linearElasticMaterial->give(Ex, gp);
            ep = CurrFt / E;
            ef = this->giveMinCrackStrainsForFullyOpenCrack(gp, ipos);
            dCoeff = ( E / ( princStress.at(ipos) / principalStrain.at(ipos) ) );
            this->giveMaterialStiffnessMatrix(Ds0, ReducedForm, SecantStiffness, gp, atTime);
            // compute reached equivalent strain
            equivStrain = this->computeCurrEquivStrain(gp, reducedStrainVector, E, atTime);
            damage = this->computeDamageCoeff(equivStrain, dCoeff, ep, ef);


            status->setDamageEpsfCoeff(ef);
            status->setDamageEpspCoeff(ep);
            status->setDamageStiffCoeff(dCoeff);
            status->setDs0Matrix(Ds0);
            status->setTempMaxEquivStrain(equivStrain);
            status->setTempDamageCoeff(damage);
            status->setTempMode(RCSDMaterialStatus :: sdMode);
        }
    } else {
        // scalar damage mode
        double ep, ef, E, dCoeff;
        double damage = 1.0;
        //int ipos;

        E = linearElasticMaterial->give(Ex, gp);
        equivStrain = this->computeCurrEquivStrain(gp, reducedStrainVector, E, atTime);
        equivStrain = max( equivStrain, status->giveTempMaxEquivStrain() );
        reducedSpaceStressVector.beProductOf(* status->giveDs0Matrix(), reducedStrainVector);
        dCoeff = status->giveDamageStiffCoeff();
        ef = status->giveDamageEpsfCoeff();
        ep = status->giveDamageEpspCoeff();
        damage = this->computeDamageCoeff(equivStrain, dCoeff, ep, ef);
        reducedSpaceStressVector.times(1.0 - damage);

        if ( form == FullForm ) {
            crossSection->giveFullCharacteristicVector(answer, gp, reducedSpaceStressVector);
        } else {
            answer = reducedSpaceStressVector;
        }

        status->letTempStressVectorBe(reducedSpaceStressVector);

        status->setTempMaxEquivStrain(equivStrain);
        status->setTempDamageCoeff(damage);
    }

    status->letTempStrainVectorBe(totalStrain);
}


void
RCSDMaterial :: giveEffectiveMaterialStiffnessMatrix(FloatMatrix &answer,
                                                     MatResponseForm form,
                                                     MatResponseMode rMode, GaussPoint *gp,
                                                     TimeStep *atTime)
//
// returns effective material stiffness matrix in full form
// for gp stress strain mode
//
{
    RCSDMaterialStatus *status = ( RCSDMaterialStatus * ) this->giveStatus(gp);

    if ( status->giveTempMode() == RCSDMaterialStatus :: rcMode ) {
        // rotating crack mode

        RCM2Material :: giveEffectiveMaterialStiffnessMatrix(answer, form, rMode, gp, atTime);
        return;
    } else {
        // rcsd mode

        if ( ( rMode == TangentStiffness ) || ( rMode == SecantStiffness ) ) {
            FloatMatrix reducedAnswer;
            IntArray mask;
            double dCoeff;

            reducedAnswer = * status->giveDs0Matrix();
            dCoeff = 1.0 - status->giveDamageCoeff();
            if ( dCoeff < RCSD_DAMAGE_EPS ) {
                dCoeff = RCSD_DAMAGE_EPS;
            }

            reducedAnswer.times(dCoeff);

            if ( form == ReducedForm ) {
                answer = reducedAnswer;
            } else {
                this->giveStressStrainMask( mask, ReducedForm, gp->giveMaterialMode() );
                answer.beSubMatrixOf(reducedAnswer, mask);
            }
        } else if ( rMode == ElasticStiffness ) {
            this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, rMode, gp, atTime);
            return;
        } else {
            _error("giveEffectiveMaterialStiffnessMatrix: usupported mode");
        }
    }
}


double
RCSDMaterial :: computeDamageCoeff(double equivStrain, double dStiffCoeff, double ep, double ef)
{
    double damage = 1. - dStiffCoeff * ( 1. - ef / equivStrain ) / ( 1. - ef / ep );
    if ( damage > 1.0 ) {
        return 1.0;
    }

    return damage;
}


double
RCSDMaterial :: computeCurrEquivStrain(GaussPoint *gp, const FloatArray &reducedTotalStrainVector, double e, TimeStep *atTime)
{
    FloatArray effStress, princEffStress, fullEffStress;
    FloatMatrix De;
    double answer = 0.0;
    int i;

    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();
    linearElasticMaterial->giveCharacteristicMatrix(De, ReducedForm, TangentStiffness, gp, atTime);
    effStress.beProductOf(De, reducedTotalStrainVector);
    crossSection->giveFullCharacteristicVector(fullEffStress, gp, effStress);

    this->computePrincipalValues(princEffStress, fullEffStress, principal_stress);
    for ( i = 1; i <= 3; i++ ) {
        answer = max( answer, macbra( princEffStress.at(i) ) );
    }

    return answer / e;
}


IRResultType
RCSDMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, SDTransitionCoeff, IFT_RCSDMaterial_sdtransitioncoeff, "sdtransitioncoeff"); // Macro

    return RCM2Material :: initializeFrom(ir);
}


double
RCSDMaterial :: give(int aProperty, GaussPoint *gp)
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
{
    if ( aProperty == pscm_SDTransitionCoeff ) {
        return this->SDTransitionCoeff;
    }

    return RCM2Material :: give(aProperty, gp);
}


int
RCSDMaterial :: checkSizeLimit(GaussPoint *gp, double charLength)
//
// checks if element size (charLength) is too big
// so that tension strength must be reduced followed
// by sudden stress drop
//
{
    double Ee, Gf, Ft, LeCrit;

    Ee = this->give(pscm_Ee, gp);
    Gf = this->give(pscm_Gf, gp);
    Ft = this->give(pscm_Ft, gp);

    LeCrit = 2.0 * Gf * Ee / ( Ft * Ft );
    return ( charLength < LeCrit );
}


double
RCSDMaterial :: computeStrength(GaussPoint *gp, double charLength)
//
// computes strength for given gp,
// which may be reduced according to length of "fracture process zone"
// to be energetically correct
//
{
    double Ee, Gf, Ft;

    Ee = this->give(pscm_Ee, gp);
    Gf = this->give(pscm_Gf, gp);
    Ft = this->give(pscm_Ft, gp);

    if ( this->checkSizeLimit(gp, charLength) ) {
        ;
    } else {
        // we reduce Ft and there is no softening but sudden drop
        Ft = sqrt(2. * Ee * Gf / charLength);
        //
        OOFEM_LOG_RELEVANT("Reducing Ft to %f in element %d, gp %d, Le %f\n",
                           Ft, gp->giveElement()->giveNumber(), gp->giveNumber(), charLength);
    }

    return Ft;
}


double
RCSDMaterial :: giveMinCrackStrainsForFullyOpenCrack(GaussPoint *gp, int i)
//
// computes MinCrackStrainsForFullyOpenCrack for given gp and i-th crack
//
{
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);
    double Le, Gf, Ft;

    Le = status->giveCharLength(i);
    Gf = this->give(pscm_Gf, gp);
    Ft = this->computeStrength(gp, Le);

    return 2.0 * Gf / ( Le * Ft );
}


/*
 * void
 * RCSDMaterial :: updateStatusForNewCrack (GaussPoint* gp, int i, double Le)
 * //
 * // updates gp status when new crack-plane i is formed with charLength Le
 * // updates Le and computes and sets minEffStrainForFullyOpenCrack
 * //
 * {
 * RCM2MaterialStatus *status = (RCM2MaterialStatus*) this -> giveStatus (gp);
 *
 * if (Le <= 0) {
 * char errMsg [80];
 * sprintf (errMsg,"Element %d returned zero char length",
 *    gp->giveElement()->giveNumber());
 * _error (errMsg);
 * }
 *
 * status -> setCharLength(i, Le);
 * status ->
 * setMinCrackStrainsForFullyOpenCrack (i, this->giveMinCrackStrainsForFullyOpenCrack(gp,i));
 * }
 */

double
RCSDMaterial :: giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp,
                                    double crackStrain, int i)
//
// returns current cracking modulus according to crackStrain for i-th
// crackplane
// now linear softening is implemented
// see also CreateStatus () function.
// softening modulus represents a relation between the normal crack strain
// rate and the normal stress rate.
//
{
    //double Ee, Gf;
    double Cf, Ft, Le, minEffStrainForFullyOpenCrack;
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);

    //
    // now we have to set proper reduced strength and softening modulus Et
    // in order to obtain energetically correct solution using the concept
    // of fracture energy Gf, which is a material constant.
    //
    //Ee = this->give(pscm_Ee);
    //Gf = this->give(pscm_Gf);
    Ft = this->give(pscm_Ft, gp);
    Le = status->giveCharLength(i);

    Ft = this->computeStrength(gp, Le);
    minEffStrainForFullyOpenCrack = this->giveMinCrackStrainsForFullyOpenCrack(gp, i);

    if ( rMode == TangentStiffness ) {
        if ( this->checkSizeLimit(gp, Le) ) {
            if ( ( crackStrain >= minEffStrainForFullyOpenCrack ) ||
                ( status->giveTempMaxCrackStrain(i) >= minEffStrainForFullyOpenCrack ) ) {
                // fully open crack - no stiffness
                Cf = 0.;
            } else if ( crackStrain >= status->giveTempMaxCrackStrain(i) ) {
                // further softening
                Cf = -Ft / minEffStrainForFullyOpenCrack;
            } else {
                // unloading or reloading regime
                Cf = Ft * ( minEffStrainForFullyOpenCrack - status->giveTempMaxCrackStrain(i) ) /
                     ( status->giveTempMaxCrackStrain(i) * minEffStrainForFullyOpenCrack );
            }
        } else {
            Cf = 0.;
        }
    } else {
        if ( this->checkSizeLimit(gp, Le) ) {
            if ( ( crackStrain >= minEffStrainForFullyOpenCrack ) ||
                ( status->giveTempMaxCrackStrain(i) >= minEffStrainForFullyOpenCrack ) ) {
                // fully open crack - no stiffness
                Cf = 0.;
            } else {
                // unloading or reloading regime
                Cf = Ft * ( minEffStrainForFullyOpenCrack - status->giveTempMaxCrackStrain(i) ) /
                     ( status->giveTempMaxCrackStrain(i) * minEffStrainForFullyOpenCrack );
            }
        } else {
            Cf = 0.;
        }
    }

    return Cf;
}



double
RCSDMaterial :: giveNormalCrackingStress(GaussPoint *gp, double crackStrain, int i)
//
// returns receivers Normal Stress in crack i  for given cracking strain
//
{
    double Cf, Ft, Le, answer, minEffStrainForFullyOpenCrack;
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);
    minEffStrainForFullyOpenCrack = this->giveMinCrackStrainsForFullyOpenCrack(gp, i);

    Cf = this->giveCrackingModulus(TangentStiffness, gp, crackStrain, i); // < 0
    Le = status->giveCharLength(i);
    Ft = this->computeStrength(gp, Le);

    if ( this->checkSizeLimit(gp, Le) ) {
        if ( ( crackStrain >= minEffStrainForFullyOpenCrack ) ||
            ( status->giveTempMaxCrackStrain(i) >= minEffStrainForFullyOpenCrack ) ) {
            // fully open crack - no stiffness
            answer = 0.;
        } else if ( crackStrain >= status->giveTempMaxCrackStrain(i) ) {
            // further softening
            answer = Ft + Cf * crackStrain;
        } else if ( crackStrain <= 0. ) {
            // crack closing
            // no stress due to cracking
            answer = 0.;
        } else {
            // unloading or reloading regime
            answer = crackStrain * Ft *
                     ( minEffStrainForFullyOpenCrack - status->giveTempMaxCrackStrain(i) ) /
                     ( status->giveTempMaxCrackStrain(i) * minEffStrainForFullyOpenCrack );
        }
    } else {
        answer = 0.;
    }

    return answer;
}




RCSDMaterialStatus :: RCSDMaterialStatus(int n, Domain *d, GaussPoint *g) :
    RCM2MaterialStatus(n, d, g), Ds0()
{
    maxEquivStrain = tempMaxEquivStrain = 0.0;
    damageCoeff = tempDamageCoeff = 1.0;
    damageStiffCoeff = depsf = depsp = 0.0;
    mode = tempMode = rcMode;
}


RCSDMaterialStatus :: ~RCSDMaterialStatus()
{ }


void
RCSDMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    int i;
    char s [ 11 ];

    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( this->giveTempMode() == rcMode ) {
        fprintf(file, "mode :rc ");

        if ( this->giveTempAlreadyCrack() ) {
            for ( i = 1; i <= 3; i++ ) {
                switch ( crackStatuses.at(i) ) {
                case pscm_NONE:
                    strcpy(s, "NONE");
                    break;
                case pscm_OPEN:
                    strcpy(s, "OPEN");
                    break;
                case pscm_SOFTENING:
                    strcpy(s, "SOFTENING");
                    break;
                case pscm_RELOADING:
                    strcpy(s, "RELOADING");
                    break;
                case pscm_UNLOADING:
                    strcpy(s, "UNLOADING");
                    break;
                default:
                    strcpy(s, "UNKNOWN");
                    break;
                }

                fprintf( file, "crack %d {status %s, normal to crackplane { %f %f %f }} ",
                        i, s, crackDirs.at(1, i), crackDirs.at(2, i), crackDirs.at(3, i) );
            }
        }
    } else {
        // sd mode
        fprintf( file, "mode :sd, damageCoeff = %f ", this->giveTempDamageCoeff() );
    }

    fprintf(file, "}\n");
}


void
RCSDMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    RCM2MaterialStatus :: initTempStatus();

    tempMaxEquivStrain = maxEquivStrain;
    tempDamageCoeff = damageCoeff;
    tempMode = mode;
}


void
RCSDMaterialStatus :: updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    RCM2MaterialStatus :: updateYourself(atTime);

    maxEquivStrain = tempMaxEquivStrain;
    damageCoeff = tempDamageCoeff;
    mode = tempMode;
}


contextIOResultType
RCSDMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = RCM2MaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream->write(& maxEquivStrain, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& damageCoeff, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& mode, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = Ds0.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
RCSDMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = RCM2MaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(& maxEquivStrain, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& damageCoeff, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& mode, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = Ds0.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK; // return succes
}

} // end namespace oofem
