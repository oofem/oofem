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

#include "rcsde.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "Materials/isolinearelasticmaterial.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Material(RCSDEMaterial);

RCSDEMaterial :: RCSDEMaterial(int n, Domain *d) : RCM2Material(n, d)
    //
    // constructor
    //
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
}


RCSDEMaterial :: ~RCSDEMaterial()
//
// destructor
//
{
    delete linearElasticMaterial;
}


void
RCSDEMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &totalStrain,
                                      TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    FloatMatrix Ds0;
    double equivStrain;
    FloatArray princStress, reducedAnswer;
    FloatArray reducedStrainVector, strainVector, principalStrain;
    FloatArray reducedSpaceStressVector;
    FloatMatrix tempCrackDirs;
    RCSDEMaterialStatus *status = static_cast< RCSDEMaterialStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedStrainVector, gp, totalStrain,
                                                tStep, VM_Total);

    StructuralMaterial :: giveFullSymVectorForm( strainVector, reducedStrainVector, gp->giveMaterialMode() );

    tempCrackDirs = status->giveTempCrackDirs();
    this->computePrincipalValDir(principalStrain, tempCrackDirs,
                                 strainVector,
                                 principal_strain);


    if ( status->giveTempMode() == RCSDEMaterialStatus :: rcMode ) {
        // rotating crack mode

        this->giveRealPrincipalStressVector3d(princStress, gp, principalStrain, tempCrackDirs, tStep);
        princStress.resize(6);
        tempCrackDirs = status->giveTempCrackDirs();
        this->transformStressVectorTo(answer, tempCrackDirs, princStress, 1);

        StructuralMaterial :: giveReducedSymVectorForm( reducedSpaceStressVector, answer, gp->giveMaterialMode() );
        status->letTempStressVectorBe(reducedSpaceStressVector);

        this->updateCrackStatus(gp, status->giveCrackStrainVector());

        StructuralMaterial :: giveReducedSymVectorForm( reducedAnswer, answer, gp->giveMaterialMode() );
        answer = reducedAnswer;

        // test if transition to scalar damage mode take place
        double minSofteningPrincStress = this->Ft, E, Le, CurrFt, Ft, Gf, Gf0, Gf1, e0, ef, ef2, damage;
        int ipos = 0;
        for ( int i = 1; i <= 3; i++ ) {
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

            Le = status->giveCharLength(ipos);
            E = linearElasticMaterial->give(Ex, gp);
            Gf = this->give(pscm_Gf, gp) / Le;
            Ft = this->computeStrength(gp, Le);
            ef = Gf / Ft;
            e0 = principalStrain.at(ipos);
            Gf0 = -CurrFt * ef * ( exp(-status->giveCrackStrain(ipos) / ef) - 1.0 ); // already disipated + 0.5*sigma0*epsilon0
            Gf1 = Gf - Gf0;

            ef2 = Gf1 / princStress.at(ipos);

            this->giveMaterialStiffnessMatrix(Ds0, SecantStiffness, gp, tStep);
            // compute reached equivalent strain
            equivStrain = this->computeCurrEquivStrain(gp, reducedStrainVector, E, tStep);
            damage = this->computeDamageCoeff(equivStrain, e0, ef2);


            status->setTransitionEpsCoeff(e0);
            status->setEpsF2Coeff(ef2);
            status->setDs0Matrix(Ds0);
            status->setTempMaxEquivStrain(equivStrain);
            status->setTempDamageCoeff(damage);
            status->setTempMode(RCSDEMaterialStatus :: sdMode);
        }
    } else {
        // scalar damage mode
        double E, e0, ef2;
        double damage;
        //int ipos;

        E = linearElasticMaterial->give(Ex, gp);
        equivStrain = this->computeCurrEquivStrain(gp, reducedStrainVector, E, tStep);
        equivStrain = max( equivStrain, status->giveTempMaxEquivStrain() );
        reducedSpaceStressVector.beProductOf(* status->giveDs0Matrix(), reducedStrainVector);
        //dCoeff = status->giveDamageStiffCoeff();
        ef2 = status->giveEpsF2Coeff();
        e0  = status->giveTransitionEpsCoeff();
        damage = this->computeDamageCoeff(equivStrain, e0, ef2);
        reducedSpaceStressVector.times(1.0 - damage);

        answer = reducedSpaceStressVector;

        status->letTempStressVectorBe(reducedSpaceStressVector);

        status->setTempMaxEquivStrain(equivStrain);
        status->setTempDamageCoeff(damage);
    }

    status->letTempStrainVectorBe(totalStrain);
}


void
RCSDEMaterial :: giveEffectiveMaterialStiffnessMatrix(FloatMatrix &answer,
                                                      MatResponseMode rMode, GaussPoint *gp,
                                                      TimeStep *tStep)
//
// returns effective material stiffness matrix in full form
// for gp stress strain mode
//
{
    RCSDEMaterialStatus *status = static_cast< RCSDEMaterialStatus * >( this->giveStatus(gp) );

    if ( status->giveTempMode() == RCSDEMaterialStatus :: rcMode ) {
        // rotating crack mode

        RCM2Material :: giveEffectiveMaterialStiffnessMatrix(answer, rMode, gp, tStep);
        return;
    } else {
        // rcsd mode

        if ( ( rMode == TangentStiffness ) || ( rMode == SecantStiffness ) ) {
            FloatMatrix reducedAnswer;
            double dCoeff;

            reducedAnswer = * status->giveDs0Matrix();
            dCoeff = 1.0 - status->giveDamageCoeff();
            if ( dCoeff < RCSDE_DAMAGE_EPS ) {
                dCoeff = RCSDE_DAMAGE_EPS;
            }

            reducedAnswer.times(dCoeff);

            answer = reducedAnswer;
        } else if ( rMode == ElasticStiffness ) {
            this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, rMode, gp, tStep);
            return;
        } else {
            OOFEM_ERROR("usupported mode");
        }
    }
}


double
RCSDEMaterial :: computeDamageCoeff(double equivStrain, double e0, double ef2)
{
    double damage = 1. - e0 / equivStrain *exp(-( equivStrain - e0 ) / ef2);
    if ( damage > 1.0 ) {
        return 1.0;
    }

    return damage;
}


double
RCSDEMaterial :: computeCurrEquivStrain(GaussPoint *gp, const FloatArray &reducedTotalStrainVector, double e, TimeStep *tStep)
{
    FloatArray effStress, princEffStress, fullEffStress;
    FloatMatrix De;
    double answer = 0.0;

    linearElasticMaterial->giveStiffnessMatrix(De, TangentStiffness, gp, tStep);
    effStress.beProductOf(De, reducedTotalStrainVector);
    StructuralMaterial :: giveFullSymVectorForm( fullEffStress, effStress, gp->giveMaterialMode() );

    this->computePrincipalValues(princEffStress, fullEffStress, principal_stress);
    for ( int i = 1; i <= 3; i++ ) {
        answer = max( answer, macbra( princEffStress.at(i) ) );
    }

    return answer / e;
}


IRResultType
RCSDEMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, SDTransitionCoeff, _IFT_RCSDEMaterial_sdtransitioncoeff);
    return RCM2Material :: initializeFrom(ir);
}


double
RCSDEMaterial :: give(int aProperty, GaussPoint *gp)
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
{
    if ( aProperty == pscm_SDTransitionCoeff ) {
        return this->SDTransitionCoeff;
    }

    return RCM2Material :: give(aProperty, gp);
}


int
RCSDEMaterial :: checkSizeLimit(GaussPoint *gp, double charLength)
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
RCSDEMaterial :: computeStrength(GaussPoint *gp, double charLength)
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

    if ( this->checkSizeLimit(gp, charLength) ) { } else {
        // we reduce Ft and there is no softening but sudden drop
        Ft = sqrt(2. * Ee * Gf / charLength);
        //
        OOFEM_LOG_RELEVANT("Reducing Ft to %f in element %d, gp %d, Le %f\n",
                           Ft, gp->giveElement()->giveNumber(), gp->giveNumber(), charLength);
    }

    return Ft;
}


double
RCSDEMaterial :: giveMinCrackStrainsForFullyOpenCrack(GaussPoint *gp, int i)
//
// computes MinCrackStrainsForFullyOpenCrack for given gp and i-th crack
//
{
    /*
     * RCM2MaterialStatus *status = (RCM2MaterialStatus*) this -> giveStatus (gp);
     * double Le, Gf, Ft;
     *
     * Le = status -> giveCharLength (i);
     * Gf = this->give(pscm_Gf);
     * Ft = this->computeStrength (gp, Le);
     *
     * return Gf/(Le*Ft); // Exponential softening
     */
    return 1.e6;
}


/*
 * void
 * RCSDEMaterial :: updateStatusForNewCrack (GaussPoint* gp, int i, double Le)
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
 * OOFEM_ERROR(errMsg);
 * }
 *
 * status -> setCharLength(i, Le);
 * status ->
 * setMinCrackStrainsForFullyOpenCrack (i, this->giveMinCrackStrainsForFullyOpenCrack(gp,i));
 * }
 */

double
RCSDEMaterial :: giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp,
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
    // double Ee, Gf;
    double Gf, Cf, Ft, Le, ef;
    RCM2MaterialStatus *status = static_cast< RCM2MaterialStatus * >( this->giveStatus(gp) );

    //
    // now we have to set proper reduced strength and softening modulus Et
    // in order to obtain energetically correct solution using the concept
    // of fracture energy Gf, which is a material constant.
    //
    //Ee = this->give(pscm_Ee);
    Gf = this->give(pscm_Gf, gp);
    Le = status->giveCharLength(i);
    Ft = this->computeStrength(gp, Le);
    //minEffStrainForFullyOpenCrack = this->giveMinCrackStrainsForFullyOpenCrack(gp,i);
    ef = Gf / ( Le * Ft );

    if ( rMode == TangentStiffness ) {
        if ( this->checkSizeLimit(gp, Le) ) {
            if ( crackStrain >= status->giveTempMaxCrackStrain(i) ) {
                // further softening
                Cf = -Ft / ef *exp(-crackStrain / ef);
            } else {
                // unloading or reloading regime
                Cf = Ft / status->giveTempMaxCrackStrain(i) * exp(-status->giveTempMaxCrackStrain(i) / ef);
            }
        } else {
            Cf = 0.;
        }
    } else {
        if ( this->checkSizeLimit(gp, Le) ) {
            Cf = Ft / status->giveTempMaxCrackStrain(i) * exp(-status->giveTempMaxCrackStrain(i) / ef);
        } else {
            Cf = 0.;
        }
    }

    return Cf;
}


double
RCSDEMaterial :: giveNormalCrackingStress(GaussPoint *gp, double crackStrain, int i)
//
// returns receivers Normal Stress in crack i  for given cracking strain
//
{
    double Gf, Ft, Le, answer, ef;
    RCM2MaterialStatus *status = static_cast< RCM2MaterialStatus * >( this->giveStatus(gp) );

    Le = status->giveCharLength(i);
    Ft = this->computeStrength(gp, Le);
    Gf = this->give(pscm_Gf, gp);
    ef = Gf / ( Le * Ft );

    if ( this->checkSizeLimit(gp, Le) ) {
        if ( crackStrain >= status->giveTempMaxCrackStrain(i) ) {
            // further softening
            answer = Ft * exp(-crackStrain / ef);
        } else {
            // crack closing
            // or unloading or reloading regime
            answer = Ft * crackStrain / status->giveTempMaxCrackStrain(i) *
            exp(-status->giveTempMaxCrackStrain(i) / ef);
        }
    } else {
        answer = 0.;
    }

    return answer;
}



RCSDEMaterialStatus :: RCSDEMaterialStatus(int n, Domain *d, GaussPoint *g) :
    RCM2MaterialStatus(n, d, g), Ds0()
{
    maxEquivStrain = tempMaxEquivStrain = 0.0;
    damageCoeff = tempDamageCoeff = 1.0;
    transitionEps = epsF2 = 0.0;
    rcsdMode = tempRcsdMode = rcMode;
}


RCSDEMaterialStatus :: ~RCSDEMaterialStatus()
{ }


void
RCSDEMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    char s [ 11 ];

    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( this->giveTempMode() == rcMode ) {
        fprintf(file, "mode :rc ");

        if ( this->giveTempAlreadyCrack() ) {
            for ( int i = 1; i <= 3; i++ ) {
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
RCSDEMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    RCM2MaterialStatus :: initTempStatus();

    tempMaxEquivStrain = maxEquivStrain;
    tempDamageCoeff = damageCoeff;
    tempRcsdMode = rcsdMode;
}


void
RCSDEMaterialStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    RCM2MaterialStatus :: updateYourself(tStep);

    maxEquivStrain = tempMaxEquivStrain;
    damageCoeff = tempDamageCoeff;
    rcsdMode = tempRcsdMode;
}


contextIOResultType
RCSDEMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
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

    if ( !stream.write(maxEquivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(damageCoeff) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    int _mode = rcsdMode;
    if ( !stream.write(_mode) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(transitionEps) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(epsF2) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = Ds0.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
RCSDEMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
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
    if ( !stream.read(maxEquivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(damageCoeff) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    int _mode;
    if ( !stream.read(_mode) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    rcsdMode = ( __rcsdModeType ) _mode;
    if ( !stream.read(transitionEps) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(epsF2) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = Ds0.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK; // return succes
}
} // end namespace oofem
