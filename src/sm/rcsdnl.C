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

#include "rcsdnl.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "nonlocalmaterialext.h"
#include "contextioerr.h"

namespace oofem {
RCSDNLMaterial :: RCSDNLMaterial(int n, Domain *d) : RCSDEMaterial(n, d), StructuralNonlocalMaterialExtensionInterface(d)
    //
    // constructor
    //
{
    //linearElasticMaterial = new IsotropicLinearElasticMaterial (n,d);
    SDTransitionCoeff2 = 0.;
    R = 0.;
}


RCSDNLMaterial :: ~RCSDNLMaterial()
//
// destructor
//
{
    //delete linearElasticMaterial;
}

Interface *
RCSDNLMaterial :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialExtensionInterfaceType ) {
        return ( StructuralNonlocalMaterialExtensionInterface * ) this;
    } else {
        return NULL;
    }
}


void
RCSDNLMaterial :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *atTime)
{
    /*  Implements the service updating local variables in given integration points,
     * which take part in nonlocal average process. Actually, no update is necessary,
     * because the value used for nonlocal averaging is strain vector used for nonlocal secant stiffness
     * computation. It is therefore necessary only to store local strain in corresponding status.
     * This service is declared at StructuralNonlocalMaterial level.
     */
    RCSDNLMaterialStatus *status = ( RCSDNLMaterialStatus * ) this->giveStatus(gp);

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);

    status->setLocalStrainVectorForAverage(strainVector);
}



void
RCSDNLMaterial :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
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
    FloatArray princStress, crackStrain, nonlocalStrain, reducedSpaceStressVector;
    FloatArray reducedNonlocStrainVector, fullNonlocStrainVector, principalStrain;
    FloatMatrix tempCrackDirs;
    RCSDNLMaterialStatus *nonlocStatus, *status = ( RCSDNLMaterialStatus * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    FloatArray nonlocalContribution;
    FloatArray reducedLocalStrainVector, localStrain;

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);
    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(atTime);

    // compute nonlocal strain increment first
    dynaList< localIntegrationRecord > *list = this->giveIPIntegrationList(gp); // !
    dynaList< localIntegrationRecord > :: iterator listIter;

    for ( listIter = list->begin(); listIter != list->end(); ++listIter ) {
        nonlocStatus = ( RCSDNLMaterialStatus * ) this->giveStatus( ( * listIter ).nearGp );
        nonlocalContribution = nonlocStatus->giveLocalStrainVectorForAverage();
        nonlocalContribution.times( ( * listIter ).weight );

        reducedNonlocStrainVector.add(nonlocalContribution);
    }

    reducedNonlocStrainVector.times( 1. / status->giveIntegrationScale() );
    this->endIPNonlocalAverage(gp);  // !

    // subtract stress independent part
    ////# this->giveStressDependentPartOfStrainVector(nonlocalStrainIncrement, gp, nonlocalTotalStrainIncrement,atTime);
    //

    reducedLocalStrainVector = totalStrain;

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(nonlocalStrain, gp, reducedNonlocStrainVector,
                                                atTime, VM_Total);
    this->giveStressDependentPartOfStrainVector(localStrain, gp, reducedLocalStrainVector,
                                                atTime, VM_Total);

    crossSection->giveFullCharacteristicVector(fullNonlocStrainVector, gp, nonlocalStrain);

    status->giveTempCrackDirs(tempCrackDirs);
    this->computePrincipalValDir(principalStrain, tempCrackDirs,
                                 fullNonlocStrainVector,
                                 principal_strain);


    if ( status->giveTempMode() == RCSDEMaterialStatus :: rcMode ) {
        // rotating crack mode

        this->giveRealPrincipalStressVector3d(princStress, gp, principalStrain, tempCrackDirs, atTime);

        ////#  //
        ////#  // this -> giveRealPrincipalStressVector3d (princStress, gp, strainIncrement, atTime);
        ////# princStress.resize (6);
        ////#  status->giveTempCrackDirs (tempCrackDirs);
        ////#  this -> transformStressVectorTo (answer, tempCrackDirs, princStress, 1);

        ////#  crossSection->giveReducedCharacteristicVector(stressIncrement, gp, answer);
        ////# stressIncrement.subtract (status -> giveStressVector());
        ////#  status -> letStressIncrementVectorBe (stressIncrement);

        status->giveCrackStrainVector(crackStrain);
        this->updateCrackStatus(gp, crackStrain);

        ////#
        this->giveMaterialStiffnessMatrix(Ds0, ReducedForm, SecantStiffness, gp, atTime);

        ////#  if (form == ReducedForm) {
        ////#   crossSection->giveReducedCharacteristicVector(reducedAnswer, gp, answer);
        ////#   answer = reducedAnswer;
        ////#  }

        // check for any currently opening crack
        int anyOpeningCrack = 0;
        for ( i = 1; i <= 3; i++ ) {
            if ( ( status->giveTempCrackStatus(i) == pscm_SOFTENING ) ||
                ( status->giveTempCrackStatus(i) == pscm_OPEN ) ) {
                anyOpeningCrack++;
            }
        }

        if ( anyOpeningCrack ) {
            // test if transition to scalar damage mode take place
            double minSofteningPrincStress = this->Ft, E, Le, CurrFt, Gf, Gf0, Gf1, e0, ef, ef2, damage;
            //  double minSofteningPrincStress = this->Ft, dCoeff, CurrFt, E, ep, ef, damage;
            int ipos = 0;
            for ( i = 1; i <= 3; i++ ) {
                if ( ( status->giveTempCrackStatus(i) == pscm_SOFTENING ) ||
                    ( status->giveTempCrackStatus(i) == pscm_OPEN ) ) {
                    if ( princStress.at(i) < minSofteningPrincStress ) {
                        minSofteningPrincStress = princStress.at(i);
                        ipos = i;
                    }
                }
            }

            CurrFt = this->computeStrength( gp, status->giveCharLength(ipos) );

            ////##
            // next pasted from rcm2:giveEffectiveMaterialStiffnessMatrix
            double G, minG, currG, princStressDis, princStrainDis;
            int ii, jj;

            minG = G = this->give(pscm_G, gp);
            for ( i = 4; i <= 6; i++ ) {
                if ( ( this->giveStressStrainComponentIndOf(FullForm, gp->giveMaterialMode(), i) ) ) {
                    if ( i == 4 ) {
                        ii = 2;
                        jj = 3;
                    } else if ( i == 5 ) {
                        ii = 1;
                        jj = 3;
                    } else if ( i == 6 ) {
                        ii = 1;
                        jj = 2;
                    } else {
                        continue;
                    }

                    princStressDis = princStress.at(ii) -
                                     princStress.at(jj);
                    princStrainDis = principalStrain.at(ii) -
                                     principalStrain.at(jj);

                    if ( fabs(princStrainDis) < rcm_SMALL_STRAIN ) {
                        currG = G;
                    } else {
                        currG = princStressDis / ( 2.0 * princStrainDis );
                    }

                    //currG = Ds0.at(indi, indi);
                    minG = min(minG, currG);
                }
            }

            ////#

            if ( ( minSofteningPrincStress <= this->SDTransitionCoeff * CurrFt ) ||
                ( minG <= this->SDTransitionCoeff2 * G ) ) {
                //   printf ("minSofteningPrincStress=%lf, CurrFt=%lf, SDTransitionCoeff=%lf",minSofteningPrincStress, CurrFt, this->SDTransitionCoeff);
                //   printf ("\nminG=%lf, G=%lf, SDTransitionCoeff2=%lf\n",minG, G, this->SDTransitionCoeff2);

                // sd transition takes place
                if ( ipos == 0 ) {
                    for ( i = 1; i <= 3; i++ ) {
                        if ( ( status->giveTempCrackStatus(i) == pscm_SOFTENING ) ||
                            ( status->giveTempCrackStatus(i) == pscm_OPEN ) ) {
                            if ( ipos == 0 ) {
                                ipos = i;
                                minSofteningPrincStress = princStress.at(i);
                            }

                            if ( princStress.at(i) < minSofteningPrincStress ) {
                                minSofteningPrincStress = princStress.at(i);
                                ipos = i;
                            }
                        }
                    }
                }

                // test for internal consistency error
                // we should switch to scalar damage, but no softening take place
                if ( ipos == 0 ) {
                    //RCSDEMaterial :: _error ("giveRealStressVector: can not switch to sd mode, while no cracking");
                    _error("giveRealStressVector: can not switch to sd mode, while no cracking");
                }

                //if (minSofteningPrincStress <= this->SDTransitionCoeff * CurrFt) printf (".");
                //else printf (":");
                //
                Le = status->giveCharLength(ipos);
                E = linearElasticMaterial->give(Ex, gp);
                Gf = this->give(pscm_Gf, gp) / Le;
                ef = this->ef;
                e0 = principalStrain.at(ipos);
                Gf0 = -CurrFt * ef * ( exp(-status->giveCrackStrain(ipos) / ef) - 1.0 ); // already disipated + 0.5*sigma0*epsilon0
                Gf1 = Gf - Gf0;

                ef2 = Gf1 / princStress.at(ipos);

                //this->giveMaterialStiffnessMatrix (Ds0, ReducedForm, SecantStiffness, gp, atTime);
                // compute reached equivalent strain
                equivStrain = this->computeCurrEquivStrain(gp, nonlocalStrain, E, atTime);
                damage = this->computeDamageCoeff(equivStrain, e0, ef2);

                //   printf ("Gf=%lf, Gf0=%lf, damage=%lf, e0=%lf, ef2=%lf, es=%lf\n",Gf, Gf0, damage,e0,ef2,equivStrain);

                status->setTransitionEpsCoeff(e0);
                status->setEpsF2Coeff(ef2);
                status->setDs0Matrix(Ds0);
                status->setTempMaxEquivStrain(equivStrain);
                status->setTempDamageCoeff(damage);
                status->setTempMode(RCSDEMaterialStatus :: sdMode);
            }
        }
    } else if ( status->giveTempMode() == RCSDEMaterialStatus :: sdMode ) {
        // scalar damage mode
        double E, e0, ef2;
        // double ep, ef, E, dCoeff;
        FloatArray reducedSpaceStressVector;
        double damage = 1.0;

        E = linearElasticMaterial->give(Ex, gp);
        equivStrain = this->computeCurrEquivStrain(gp, nonlocalStrain, E, atTime);
        equivStrain = max( equivStrain, status->giveTempMaxEquivStrain() );
        ////# reducedSpaceStressVector.beProductOf (*status->giveDs0Matrix(), reducedNonlocStrainVector);
        ef2 = status->giveEpsF2Coeff();
        e0  = status->giveTransitionEpsCoeff();
        damage = this->computeDamageCoeff(equivStrain, e0, ef2);

        //  dCoeff = status->giveDamageStiffCoeff ();
        //  ef = status->giveDamageEpsfCoeff();
        //  ep = status->giveDamageEpspCoeff();
        //  damage = this->computeDamageCoeff (equivStrain, dCoeff, ep, ef);
        ////#
        Ds0 = * status->giveDs0Matrix();
        Ds0.times(1.0 - damage);
        ////#
        ////# reducedSpaceStressVector.times (1.0 - damage);

        ////#  if (form == FullForm) {
        ////#   crossSection->giveFullCharacteristicVector(answer, gp, reducedSpaceStressVector);
        ////#  } else {
        ////#   answer = reducedSpaceStressVector;
        ////#  }

        ////#  stressIncrement = reducedSpaceStressVector;
        ////#  stressIncrement.subtract (status -> giveStressVector());
        ////#  status -> letStressIncrementVectorBe (stressIncrement);

        status->setTempMaxEquivStrain(equivStrain);
        status->setTempDamageCoeff(damage);
    }

    ////# common part
    reducedSpaceStressVector.beProductOf(Ds0, localStrain);

    if ( form == FullForm ) {
        crossSection->giveFullCharacteristicVector(answer, gp, reducedSpaceStressVector);
    } else {
        answer = reducedSpaceStressVector;
    }

    status->letTempStressVectorBe(reducedSpaceStressVector);

    status->letTempStrainVectorBe(totalStrain);
    status->setTempNonlocalStrainVector(reducedNonlocStrainVector);
}


IRResultType
RCSDNLMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    //RCSDEMaterial::instanciateFrom (ir);
    this->giveLinearElasticMaterial()->initializeFrom(ir);
    IR_GIVE_FIELD(ir, Ft, IFT_RCSDNLMaterial_ft, "ft"); // Macro
    StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, SDTransitionCoeff, IFT_RCSDNLMaterial_sdtransitioncoeff, "sdtransitioncoeff"); // Macro
    IR_GIVE_FIELD(ir, SDTransitionCoeff2, IFT_RCSDNLMaterial_sdtransitioncoeff2, "sdtransitioncoeff2"); // Macro
    if ( SDTransitionCoeff2 > 1.0 ) {
        SDTransitionCoeff2 = 1.0;
    }

    if ( SDTransitionCoeff2 < 0.0 ) {
        SDTransitionCoeff2 = 0.0;
    }

    IR_GIVE_FIELD(ir, R, IFT_RCSDNLMaterial_r, "r"); // Macro
    if ( R < 0.0 ) {
        R = 0.0;
    }

    if ( ir->hasField(IFT_RCSDNLMaterial_ef, "ef") ) { // if ef is specified, Gf is computed acordingly
        IR_GIVE_FIELD(ir, this->ef, IFT_RCSDNLMaterial_ef, "ef"); // Macro
        this->Gf = this->Ft * this->ef;
    } else if ( ir->hasField(IFT_RCSDNLMaterial_gf, "gf") ) { // otherwise if Gf is specified, ef is computed acordingly
        IR_GIVE_FIELD(ir, this->Gf, IFT_RCSDNLMaterial_gf, "gf"); // Macro
        this->ef = this->Gf / this->Ft;
    } else {
        _error("initializeFrom: cannont determine Gf and ef from input data");
    }

    return IRRT_OK;
}


double
RCSDNLMaterial :: giveMinCrackStrainsForFullyOpenCrack(GaussPoint *gp, int i)
{
    return 1.e6; //this->ef;
}


double
RCSDNLMaterial :: computeWeightFunction(const FloatArray &src, const FloatArray &coord)
{
    // Bell shaped function decaying with the distance.

    double dist = src.distance(coord);

    if ( ( dist >= 0. ) && ( dist <= this->R ) ) {
        double help = ( 1. - dist * dist / ( R * R ) );
        return help * help;
    }

    return 0.0;
}
/*
 * void
 * RCSDNLMaterial :: updateStatusForNewCrack (GaussPoint* gp, int i, double Le)
 * //
 * // updates gp status when new crack-plane i is formed with charLength Le
 * // updates Le and computes and sets minEffStrainForFullyOpenCrack
 * //
 * {
 * RCSDNLMaterialStatus *status = (RCSDNLMaterialStatus*) this -> giveStatus (gp);
 *
 * if (Le <= 0) {
 * char errMsg [80];
 * sprintf (errMsg,"Element %d returned zero char length",
 *    gp->giveElement()->giveNumber());
 * RCSDMaterial::_error (errMsg);
 * }
 *
 * //status -> setCharLength(i, Le);
 * status -> setCharLength(i, 1.0);
 * status ->
 * setMinCrackStrainsForFullyOpenCrack (i, this->giveMinCrackStrainsForFullyOpenCrack(gp,i));
 * }
 */



RCSDNLMaterialStatus :: RCSDNLMaterialStatus(int n, Domain *d, GaussPoint *g) :
    RCSDEMaterialStatus(n, d, g), StructuralNonlocalMaterialStatusExtensionInterface(), nonlocalStrainVector(),
    tempNonlocalStrainVector(), localStrainVectorForAverage()
{
    nonlocalStrainVector.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->
                                giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );

    localStrainVectorForAverage.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->
                                       giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
}


RCSDNLMaterialStatus :: ~RCSDNLMaterialStatus()
{ }


void
RCSDNLMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    int i, n;
    FloatArray helpVec;

    RCSDEMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "nonlocstatus { ");
    fprintf(file, "  nonloc strains ");
    ( ( StructuralCrossSection * )
     gp->giveCrossSection() )->giveFullCharacteristicVector(helpVec, gp, nonlocalStrainVector);
    n = helpVec.giveSize();
    for ( i = 1; i <= n; i++ ) {
        fprintf( file, " % .4e", helpVec.at(i) );
    }

    fprintf(file, "}\n");
}


void
RCSDNLMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    RCSDEMaterialStatus :: initTempStatus();

    if ( nonlocalStrainVector.giveSize() == 0 ) {
        nonlocalStrainVector.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->
                                    giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
    }

    if ( localStrainVectorForAverage.giveSize() == 0 ) {
        localStrainVectorForAverage.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->
                                           giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
    }

    tempNonlocalStrainVector = nonlocalStrainVector;
}


void
RCSDNLMaterialStatus :: updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    RCSDEMaterialStatus :: updateYourself(atTime);
    nonlocalStrainVector = tempNonlocalStrainVector;
}


contextIOResultType
RCSDNLMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = RCSDEMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( ( iores = nonlocalStrainVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
RCSDNLMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = RCSDEMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data

    if ( ( iores = nonlocalStrainVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK; // return success
}


Interface *
RCSDNLMaterialStatus :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return this;
    } else {
        return NULL;
    }
}


#ifdef __PARALLEL_MODE
int
RCSDNLMaterial :: packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
{
    RCSDNLMaterialStatus *status = ( RCSDNLMaterialStatus * ) this->giveStatus(ip);

    this->buildNonlocalPointTable(ip);
    this->updateDomainBeforeNonlocAverage(stepN);

    return status->giveLocalStrainVectorForAverage().packToCommBuffer(buff);
}

int
RCSDNLMaterial :: unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
{
    int result;
    RCSDNLMaterialStatus *status = ( RCSDNLMaterialStatus * ) this->giveStatus(ip);
    FloatArray localStrainVectorForAverage;

    result = localStrainVectorForAverage.unpackFromCommBuffer(buff);
    status->setLocalStrainVectorForAverage(localStrainVectorForAverage);
    return result;
}

int
RCSDNLMaterial :: estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip)
{
    //
    // Note: status localStrainVectorForAverage memeber must be properly sized!
    //
    RCSDNLMaterialStatus *status = ( RCSDNLMaterialStatus * ) this->giveStatus(ip);

    return status->giveLocalStrainVectorForAverage().givePackSize(buff);
}

#endif
} // end namespace oofem
