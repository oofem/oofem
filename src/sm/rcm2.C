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

#include "rcm2.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "contextioerr.h"
#include "mathfem.h"

#ifndef __MAKEDEPEND
 #include <cstring>
#endif

namespace oofem {
#define rcm_STRESSRELERROR 1.e-5
#define rcm_RESIDUALSTIFFFACTOR 1.e-3

RCM2Material :: RCM2Material(int n, Domain *d) : StructuralMaterial(n, d)
//
// constructor
//
{
    linearElasticMaterial = NULL;
    Gf = Ft = 0.;
}


RCM2Material :: ~RCM2Material()
//
// destructor
//
{
    //delete linearElasticMaterial;
}

int
RCM2Material :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    if ( ( mode == _3dMat ) || ( mode == _PlaneStress ) ||
        ( mode == _PlaneStrain ) || ( mode == _1dMat ) ||
        ( mode == _2dPlateLayer ) || ( mode == _2dBeamLayer ) ||
        ( mode == _3dShellLayer ) ) {
        return 1;
    }

    return 0;
}


/*
 * void
 * RCM2Material :: computeTrialStressIncrement (FloatArray& answer, GaussPoint *gp,
 *                     const FloatArray &strainIncrement,
 *                          TimeStep* atTime)
 * //
 * // returns trial stress according to strainIncrement
 * //
 * {
 * FloatMatrix def;
 *
 * this -> giveEffectiveMaterialStiffnessMatrix (def, TangentStiffness,gp,atTime);
 * answer.beProductOf (def, strainIncrement);
 * }
 */


void
RCM2Material :: giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseForm form,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *atTime)
//
// computes full constitutive matrix for case of gp stress-strain state.
//
{
    /*
     * FloatMatrix *def;
     * RCM2MaterialStatus *status =
     * (RCM2MaterialStatus*) this -> giveStatus (gp);
     *
     * CrossSection *crossSection = gp -> giveElement()->giveCrossSection();
     */
    this->giveEffectiveMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
    // def is full matrix for current gp stress strain mode.
}


void
RCM2Material :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                     const FloatArray &totalStrain,
                                     TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    FloatArray princStress, reducedStressVector, crackStrain, reducedAnswer;
    FloatArray reducedTotalStrainVector, principalStrain, strainVector;
    FloatMatrix tempCrackDirs;
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, atTime, VM_Total);
    //

    crossSection->giveFullCharacteristicVector(strainVector, gp, reducedTotalStrainVector);

    status->giveTempCrackDirs(tempCrackDirs);
    this->computePrincipalValDir(principalStrain, tempCrackDirs,
                                 strainVector,
                                 principal_strain);

    this->giveRealPrincipalStressVector3d(princStress, gp, principalStrain, tempCrackDirs, atTime);

    princStress.resize(6);
    status->giveTempCrackDirs(tempCrackDirs);
    this->transformStressVectorTo(answer, tempCrackDirs, princStress, 1);

    status->letTempStrainVectorBe(totalStrain);
    crossSection->giveReducedCharacteristicVector(reducedStressVector, gp, answer);
    status->letTempStressVectorBe(reducedStressVector);
    status->giveCrackStrainVector(crackStrain);
    this->updateCrackStatus(gp, crackStrain);

    if ( form == FullForm ) {
        return;
    }

    crossSection->giveReducedCharacteristicVector(reducedAnswer, gp, answer);
    answer = reducedAnswer;
}


void
RCM2Material ::  giveRealPrincipalStressVector3d(FloatArray &answer, GaussPoint *gp,
                                                 FloatArray &principalStrain,
                                                 FloatMatrix &tempCrackDirs,
                                                 TimeStep *atTime)
//
// returns real principal stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
// updates principal strain and stress of the receiver's status.
//
{
    int i, iter, ind;
    double maxErr;
    FloatArray crackStrainVector, reducedTotalStrainVector;
    FloatArray strainIncrement, crackStrainIterativeIncrement;
    FloatArray prevPrincipalStrain;
    FloatArray dSigma;
    FloatArray elastStrain, sigmaEl, sigmaCr(3);
    FloatArray fullDSigma;
    IntArray activatedCracks, crackMapping;
    FloatMatrix dcr, de, decr, fullDecr, crackDirs;
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);

    /*
     * if (status -> giveStressVector() == NULL) status->letStressVectorBe(new FloatArray(this->giveSizeOfReducedStressStrainVector(gp->giveMaterialMode())));
     * if (status -> giveStrainVector() == NULL) status->letStrainVectorBe(new FloatArray(this->giveSizeOfReducedStressStrainVector(gp->giveMaterialMode())));
     * // if (status -> givePlasticStrainVector() == NULL) status->letPlasticStrainVectorBe(new FloatArray(6));
     * if (status -> giveStressIncrementVector() == NULL) status->letStressIncrementVectorBe(new FloatArray(this->giveSizeOfReducedStressStrainVector(gp->giveMaterialMode())));
     * if (status -> giveStrainIncrementVector() == NULL) status->letStrainIncrementVectorBe(new FloatArray(this->giveSizeOfReducedStressStrainVector(gp->giveMaterialMode())));
     * // if (status -> givePlasticStrainIncrementVector() == NULL) status->letPlasticStrainIncrementVectorBe(new FloatArray(6));
     */

    /*
     * // totalStressVector = gp -> giveStressVector()->GiveCopy();
     * reducedTotalStrainVector = status -> giveStrainVector();
     * reducedTotalStrainVector.add(fullStrainIncrement);
     * crossSection->giveFullCharacteristicVector(totalStrainVector, gp, reducedTotalStrainVector);
     * //delete reducedTotalStrainVector;
     * // plasticStrainVector = status -> givePlasticStrainVector()->GiveCopy();
     *
     *
     * // already cracked - next directions are determined
     * // according to principal strain directions
     * status->giveTempCrackDirs(tempCrackDirs);
     * this->computePrincipalValDir (principalStrain, tempCrackDirs,
     *              totalStrainVector,
     *              principal_strain);
     * status->letTempCrackDirsBe (tempCrackDirs);
     */
    status->giveCrackStrainVector(crackStrainVector); // local one
    status->giveCrackDirs(crackDirs);
    if ( principalStrain.containsOnlyZeroes() ) {
        // keep old principal values
        status->letTempCrackDirsBe(crackDirs);
    } else {
        this->sortPrincDirAndValCloseTo(& principalStrain,
                                        & tempCrackDirs, & crackDirs);
        status->letTempCrackDirsBe(tempCrackDirs);
    }

    // compute de in local system
    // for iso materials no transformation if stiffness required
    //
    // local strain increment
    status->givePrevPrincStrainVector(prevPrincipalStrain);
    strainIncrement.beDifferenceOf(principalStrain, prevPrincipalStrain);
    status->letPrincipalStrainVectorBe(principalStrain);

    this->giveNormalElasticStiffnessMatrix(de, FullForm, TangentStiffness,
                                           gp, atTime, tempCrackDirs);
    //
    // construct mapping matrix of active cracks
    // this mapping will dynamically change as
    // some crack can unlo or reload
    //
    this->updateActiveCrackMap(gp);
    status->giveCrackMap(crackMapping);
    // start iteration until stress computed from elastic increment
    // is equal to stress computed from cracking strain increment
    // we do this computation in reduced stress strain space
    dSigma.resize(0);
    for ( iter = 1; iter <= 20; iter++ ) {
        //
        // first check if already cracked
        //
        if ( status->giveNumberOfTempActiveCracks() ) {
            // active crack exist
            this->giveCrackedStiffnessMatrix(dcr, TangentStiffness, gp, atTime);
            fullDecr = de;
            fullDecr.add(dcr);
            decr.beSubMatrixOf(fullDecr, crackMapping);

            if ( dSigma.giveSize() == 0 ) {
                fullDSigma.beProductOf(de, strainIncrement);
                dSigma.beSubArrayOf(fullDSigma, crackMapping);
            }

            decr.solveForRhs(dSigma, crackStrainIterativeIncrement);
            for ( i = 1; i <= 3; i++ ) {
                if ( ( ind = crackMapping.at(i) ) ) {
                    crackStrainVector.at(i) += crackStrainIterativeIncrement.at(ind);
                }
            }

            // check for crack closing, updates also cracking map
            this->checkIfClosedCracks(gp, crackStrainVector, crackMapping);

            // elastic strain component
            elastStrain.beDifferenceOf(principalStrain, crackStrainVector);
            sigmaEl.beProductOf(de, elastStrain);

            // Stress in cracks
            for ( i = 1; i <= 3; i++ ) {
                if ( crackMapping.at(i) ) {
                    sigmaCr.at(i) = giveNormalCrackingStress(gp, crackStrainVector.at(i), i);
                }
            }

            // update status
            status->letCrackStrainVectorBe(crackStrainVector);
        } else {
            //
            // no active crack exist - elastic behaviour
            //
            elastStrain.beDifferenceOf(principalStrain, crackStrainVector);
            sigmaEl.beProductOf(de, elastStrain);
            sigmaCr.zero();
        }

        // check for new cracks
        // and update crack map if necessary
        // when we update map, we need to add new crack at end
        // because sigmaCr is build
        this->checkForNewActiveCracks(activatedCracks, gp, crackStrainVector,
                                      sigmaEl, sigmaCr, principalStrain);
        if ( activatedCracks.giveSize() ) {
            // update crack map also
            this->updateActiveCrackMap(gp, & activatedCracks);
            status->giveCrackMap(crackMapping); // update crackMap
        }

        //

        // compute unbalanced stress
        // dSigma = sigmaEl - sigmaCr for active cracks
        fullDSigma = sigmaEl;
        fullDSigma.subtract(sigmaCr);
        dSigma.beSubArrayOf(fullDSigma, crackMapping);
        // find max error in dSigma
        // if max err < allovedErr -> stop iteration
        // allowed Err is computed relative to Ft;

        // check only for active cracks
        maxErr = 0.;
        for ( i = 1; i <= dSigma.giveSize(); i++ ) {
            if ( fabs( dSigma.at(i) ) > maxErr ) {
                maxErr = fabs( dSigma.at(i) );
            }
        }

        if ( maxErr < rcm_STRESSRELERROR * this->give(pscm_Ft, gp) ) {
            status->letPrincipalStressVectorBe(sigmaEl);
            answer = sigmaEl;
            return;
        }
    } // loop

    // convergence not reached
    _error("GiveRealStressVector3d - convergence not reached");
}


void
RCM2Material :: initTempStatus(GaussPoint *gp)
//
// Initialize MatStatus (respective it's temporary variables at the begining
// of integrating incremental constitutive relations) to correct values
//
{
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);

    status->initTempStatus();
}


void
RCM2Material :: checkForNewActiveCracks(IntArray &answer, GaussPoint *gp,
                                        const FloatArray &crackStrain,
                                        const FloatArray &princStressVector,
                                        FloatArray &crackStressVector,
                                        const FloatArray &princStrainVector)
//
// returns int_array flag showing if some crack
// is newly activated or
// closed crack is reopened
// return 0 if no crack is activated or reactivated.
// modifies crackStressVector for newly activated crack.
//
{
    double initStress, Le = 0.0;
    int i, j, upd, activationFlag = 0;
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);
    FloatArray localStress;
    FloatMatrix tempCrackDirs;
    IntArray crackMap;

    answer.resize(3);
    answer.zero();
    status->giveCrackMap(crackMap);
    localStress = princStressVector;
    //
    // local stress is updated according to reached local crack strain
    //
    for ( i = 1; i <= 3; i++ ) { // loop over each possible crack plane
        // test previous status of each possible crack plane
        upd = 0;
        if ( ( crackMap.at(i) == 0 ) &&
            ( this->giveStressStrainComponentIndOf(FullForm, gp->giveMaterialMode(), i) ) ) {
            if ( status->giveTempMaxCrackStrain(i) > 0. ) {
                // if (status->giveTempCrackStatus()->at(i) != pscm_NONE) {
                //
                // previously cracked direction
                //
                initStress = 0.;
            } else {
                // newer cracked, so we compute principal stresses
                // and compare them with reduced tension strength
                // as a criterion for crack initiation

                FloatArray crackPlaneNormal(3);
                status->giveTempCrackDirs(tempCrackDirs);
                for ( j = 1; j <= 3; j++ ) {
                    crackPlaneNormal.at(j) = tempCrackDirs.at(j, i);
                }

                // Le = gp->giveElement()->giveCharacteristicLenght (gp, &crackPlaneNormal);
                Le = this->giveCharacteristicElementLenght(gp, crackPlaneNormal);
                initStress = this->computeStrength(gp, Le);
                upd = 1;
            }

            if ( localStress.at(i) > initStress ) {
                crackStressVector.at(i) = initStress;
                answer.at(i) = 1;
                activationFlag = 1;
                if ( upd ) {
                    this->updateStatusForNewCrack(gp, i, Le);
                }
            }
        } // end of tested crack

    } // end of loop over are possible directions

    if ( activationFlag ) {
        return;
    }

    answer.resize(0);
}


double
RCM2Material :: giveCharacteristicElementLenght(GaussPoint *gp, const FloatArray &crackPlaneNormal) {
    // returns characteristic element length for given crack plane normal
    return gp->giveElement()->giveCharacteristicLenght(gp, crackPlaneNormal);
}


void
RCM2Material :: updateStatusForNewCrack(GaussPoint *gp, int i, double Le)
//
// updates gp status when new crack-plane i is formed with charLength Le
// updates Le and computes and sets minEffStrainForFullyOpenCrack
//
{
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);

    if ( Le <= 0 ) {
        _error2( "Element %d returned zero char length", gp->giveElement()->giveNumber() );
    }

    status->setCharLength(i, Le);
    //status ->  setMinCrackStrainsForFullyOpenCrack (i, this->giveMinCrackStrainsForFullyOpenCrack(gp,i));
}


void
RCM2Material :: updateCrackStatus(GaussPoint *gp, const FloatArray &crackStrain)
//
// updates gp records and MatStatus due to cracking.
//
// Updates MatStatus (respective it's temporary variables) to current
// reached status during integrating incremental constitutive relations
// Temporary variables are used, because we may integrate constitutive
// realtions many times for different strainIncrement in order to
// reach equilibrium state, so we don't want to change other variables
// which describing previously reached equlibrium. After a new equilibrium
// is reached, this->updateYourself is called which invokes matStatus->
// updateYourself(), which copies temporary variables to variables
// describing equilibrium.
//
//
{
    int i;
    double minCrackStrainsForFullyOpenCrack;
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);
    IntArray crackMap;

    status->giveCrackMap(crackMap);

    // check if material previously cracked
    // and compute possible crack planes
    // or if newer cracked, so we compute principal stresses
    // (for iso mat. coincide with princ strains)
    // and compare them with reduced tension strength
    // as a criterion for crack initiation

    // principalStrain = status->givePrincipalStrainVector();
    // transform stresses to principal directions of strains
    //
    // local stress is updated according to reached local crack strain
    //
    for ( i = 1; i <= 3; i++ ) { // loop over each possible crack plane
        if ( ( crackMap.at(i) != 0 ) &&
            ( this->giveStressStrainComponentIndOf(FullForm, gp->giveMaterialMode(), i) ) ) {
            if ( status->giveTempMaxCrackStrain(i) < crackStrain.at(i) ) {
                status->setTempMaxCrackStrain( i,  crackStrain.at(i) );
            }

            minCrackStrainsForFullyOpenCrack = this->giveMinCrackStrainsForFullyOpenCrack(gp, i);
            //if ((crackStrain.at(i) >= status->giveMinCrackStrainsForFullyOpenCrack(i)) &&
            if ( ( crackStrain.at(i) >= minCrackStrainsForFullyOpenCrack ) &&
                ( crackStrain.at(i) >= status->giveTempMaxCrackStrain(i) ) ) {
                //
                // fully open crack
                //
                status->setTempCrackStatus(i, pscm_OPEN);
            } else if ( crackStrain.at(i) >= status->giveTempMaxCrackStrain(i) ) {
                //
                // further softening of crack
                //
                status->setTempCrackStatus(i, pscm_SOFTENING);
                // status->giveTempReachedSofteningStress()->at(i) = localStress->at(i);
            } else if ( crackStrain.at(i) <= 0. ) {
                if ( status->giveTempCrackStatus(i) != pscm_NONE ) {
                    // previously active crack becomes closed
                    status->setTempCrackStatus(i, pscm_CLOSED);
                }
            } else {
                //
                // crack unloading or reloading
                //
                status->setTempCrackStatus(i, pscm_UNLOADING);
            }
        } // end for possible direction

    } // end loop over prin directions

}


void
RCM2Material :: checkIfClosedCracks(GaussPoint *gp, FloatArray &crackStrainVector,
                                    IntArray &crackMap)
//
// Check if crack closing occurs
// if yes updates crackStrainVector and gp-status accordingly
//
{
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);
    int i, isClosing = 0;

    for ( i = 1; i <= 3; i++ ) {
        if ( crackMap.at(i) ) {
            if ( crackStrainVector.at(i) < 0. ) {
                // crack closing occur
                crackStrainVector.at(i) = 0.;
                crackMap.at(i) = 0;
                // status->giveTempCrackStatus()->at(i) = pscm_CLOSED;
                isClosing = 1;
            }
        }
    }

    status->letCrackMapBe(crackMap);
    if ( isClosing ) {
        this->updateActiveCrackMap(gp);
    }

    status->giveCrackMap(crackMap); // update crack Map
}


void
RCM2Material :: giveNormalElasticStiffnessMatrix(FloatMatrix &answer,
                                                 MatResponseForm form, MatResponseMode rMode,
                                                 GaussPoint *gp, TimeStep *atTime,
                                                 const FloatMatrix &dir)
//
// return Elastic Stiffness matrix for normal Stresses
// taking into account the directions of normal stresses
// (not supported now)
//
{
    StructuralMaterial *lMat = static_cast< StructuralMaterial * >( this->giveLinearElasticMaterial() );
    FloatMatrix de, fullAnswer(3, 3);
    IntArray mask;
    int i, j, sd;

    lMat->giveCharacteristicMatrix(de, FullForm, rMode, gp, atTime);
    // copy first 3x3 submatrix to answer
    for ( i = 1; i <= 3; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            fullAnswer.at(i, j) = de.at(i, j);
        }
    }

    if ( form == FullForm ) { // 3x3 full form required
        answer = fullAnswer;
    } else {
        // reduced form for only
        // nonzero normal stresses
        //
        // first find spatial dimension of problem
        sd = 0;
        for ( i = 1; i <= 3; i++ ) {
            if ( ( this->giveStressStrainComponentIndOf(FullForm, gp->giveMaterialMode(), i) ) ) {
                sd++;
            }
        }

        answer.resize(sd, sd);
        answer.zero();

        // copy fullAnswer to reduced one
        this->giveStressStrainMask( mask, FullForm, gp->giveMaterialMode() );
        for ( i = 1; i <= 3; i++ ) {
            for ( j = 1; j <= 3; j++ ) {
                if ( mask.at(i) && mask.at(j) ) {
                    answer.at( mask.at(i), mask.at(j) ) = fullAnswer.at(i, j);
                }
            }
        }
        return;
    }
}


void
RCM2Material :: giveEffectiveMaterialStiffnessMatrix(FloatMatrix &answer,
                                                     MatResponseForm form,
                                                     MatResponseMode rMode, GaussPoint *gp,
                                                     TimeStep *atTime)
//
// returns effective material stiffness matrix in full form
// for gp stress strain mode
//
{
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);
    StructuralMaterial *lMat = static_cast< StructuralMaterial * >( this->giveLinearElasticMaterial() );
    int numberOfActiveCracks = status->giveNumberOfTempActiveCracks();
    int i, j, indi, indj, ii, jj;
    double G, princStressDis, princStrainDis;
    FloatMatrix de, invDe, compliance, dcr, d, df, t, tt, tempCrackDirs;
    FloatArray principalStressVector, principalStrainVector;
    IntArray mask;

    if ( ( rMode == ElasticStiffness ) || ( numberOfActiveCracks == 0 ) ) {
        lMat->giveCharacteristicMatrix(answer, form, rMode, gp, atTime);
        return;
    }

    // this->updateActiveCrackMap(gp) must be done after restart.
    this->updateActiveCrackMap(gp);
    status->giveTempCrackDirs(tempCrackDirs);
    this->giveNormalElasticStiffnessMatrix(de, ReducedForm, rMode, gp, atTime,
                                           tempCrackDirs);
    invDe.beInverseOf(de);
    this->giveCrackedStiffnessMatrix(dcr, rMode, gp, atTime);
    this->giveStressStrainMask( mask, ReducedForm, gp->giveMaterialMode() );
    compliance.resize( mask.giveSize(), mask.giveSize() );

    // we will set
    // first we set compliances for normal streses in
    // local coordinate system defined by crackplane
    for ( i = 1; i <= 3; i++ ) {
        if ( ( indi = this->giveStressStrainComponentIndOf(FullForm, gp->giveMaterialMode(), i) ) ) {
            for ( j = 1; j <= 3; j++ ) {
                if ( ( indj = this->giveStressStrainComponentIndOf(FullForm, gp->giveMaterialMode(), j) ) ) {
                    compliance.at(indi, indj) += invDe.at(i, j);
                }
            }

            if ( status->isCrackActive(i) ) {
                if ( dcr.at(i, i) <= 1.e-8 ) {
                    compliance.at(indi, indi) *= rcm2_BIGNUMBER;
                } else {
                    compliance.at(indi, indi) += 1. / dcr.at(i, i);
                }
            }
        }
    }

    status->getPrincipalStressVector(principalStressVector);
    status->getPrincipalStrainVector(principalStrainVector);

    // now remain to set shears
    G = this->give(pscm_G, gp);
    for ( i = 4; i <= 6; i++ ) {
        if ( ( indi = this->giveStressStrainComponentIndOf(FullForm, gp->giveMaterialMode(), i) ) ) {
            if ( i == 4 ) {
                ii = 2;
                jj = 3;
            } else if ( i == 5 ) {
                ii = 1;
                jj = 3;
            } else {
                ii = 1;
                jj = 2;
            }

            princStressDis = principalStressVector.at(ii) -
                             principalStressVector.at(jj);
            princStrainDis = principalStrainVector.at(ii) -
                             principalStrainVector.at(jj);
            if ( fabs(princStrainDis) < rcm_SMALL_STRAIN ) {
                compliance.at(indi, indi) = 1. / G;
            } else if ( fabs(princStressDis) < 1.e-8 ) {
                compliance.at(indi, indi) = rcm2_BIGNUMBER;
            } else {
                compliance.at(indi, indi) = 2 * princStrainDis / princStressDis;
            }
        }
    }

    // now we invert compliance to get stiffness in reduced space
    d.beInverseOf(compliance);
    // delete compliance;
    //
    // now let d to grow to Full Format
    //
    this->giveStressStrainMask( mask, ReducedForm, gp->giveMaterialMode() );
    df.beSubMatrixOfSizeOf(d, mask, 6);
    //
    // final step - transform stiffnes to global c.s
    //

    this->giveStressVectorTranformationMtrx(t, tempCrackDirs, 1);
    tt.beTranspositionOf(t);
    df.rotatedWith(tt);


    if ( form == FullForm ) {
        answer = df;
    } else { // reduced form asked
        this->giveStressStrainMask( mask, FullForm, gp->giveMaterialMode() );
        answer.beSubMatrixOf(df, mask);
    }
}


void
RCM2Material :: giveCrackedStiffnessMatrix(FloatMatrix &answer,
                                           MatResponseMode rMode,
                                           GaussPoint *gp,
                                           TimeStep *atTime)
//
//
// Returns material incremental stiffness matrix for cracked concrete.
// This matrix is composed from submatrix, each of them corresponds to
// one active crack in material point.
// when constructing submatrix, following assumptions are made:
//
// A salient characteriastic of crack formation concerns the fact that in most general case
// of tree-dimensional solid only 3 out of 6 components of the crack strain rate vector
// are possibly non-zero.(the normal and two shear strain rates).
// We therefore assume that the stress-strain law for the crack has a structure
// such that the other strains rate components vanish. Moreover we assume that
// the novanishing strain rate components are only related to corresponding
// stress rate components (submatrix has dimensions 3x3).
//
// if strainIncrement is defined (not null) then we take care about possible unlo&reloading
// we don't teke care about possible cracking (or non-linear softening) during strain increment
// this correction is made in this -> updateCrackStatus  (gp);
{
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);
    int i;
    int numberOfActiveCracks = status->giveNumberOfTempActiveCracks();
    IntArray crackMap;

    status->giveCrackMap(crackMap);
    if ( numberOfActiveCracks == 0 ) {
        answer.resize(0, 0);
        return;
    }

    answer.resize(3, 3);
    answer.zero();

    // loop over each active crack plane
    for ( i = 1; i <= 3; i++ ) {
        if ( crackMap.at(i) ) {
            // obtain incremental law for one crack
            answer.at(i, i) = this->giveCrackingModulus(rMode, gp,
                                                        status->giveCrackStrain(i),
                                                        i);
        }
    }
}


void
RCM2Material :: updateActiveCrackMap(GaussPoint *gp, const IntArray *activatedCracks)
//
//
// updates mapping matrix of active cracks
//
{
    int i, indx = 1;
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);
    IntArray crackMap;
    status->giveCrackMap(crackMap);

    //if (crackMap == NULL) _error ("updateActiveCrackMap: NULL pointer encountered");

    for ( i = 1; i <= 3; i++ ) {
        if ( status->isCrackActive(i) ) {
            crackMap.at(i) = indx++;
        } else if ( activatedCracks ) {
            if ( ( activatedCracks->at(i) != 0 ) ) {
                crackMap.at(i) = indx++;
            }
        } else {
            crackMap.at(i) = 0;
        }
    }

    // store modified map into status
    status->letCrackMapBe(crackMap);
}


IRResultType
RCM2Material :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, Gf, IFT_RCM2Material_gf, "gf"); // Macro
    IR_GIVE_FIELD(ir, Ft, IFT_RCM2Material_ft, "ft"); // Macro

    this->giveLinearElasticMaterial()->initializeFrom(ir);

    return IRRT_OK;
}

double
RCM2Material :: give(int aProperty, GaussPoint *gp)
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
{
    double value = 0.0;

    if ( aProperty == pscm_Gf ) {
        return this->Gf;
    }

    // if (aProperty == pscm_Beta) return this->beta;
    if ( aProperty == pscm_Ft ) {
        return this->Ft;
    }

    if ( aProperty == pscm_Ee ) {
        aProperty = 'E';
    }

    if ( aProperty == pscm_G ) {
        aProperty = 'G';
    }

    if ( propertyDictionary->includes(aProperty) ) {
        value = propertyDictionary->at(aProperty);
    } else {
        if ( linearElasticMaterial ) {
            value = this->linearElasticMaterial->give(aProperty, gp);
        } else {
            _error("give: property not defined");
        }
    }

    return value;
}


int
RCM2Material :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(aGaussPoint);
    if ( type == IST_CrackedFlag ) {
        answer.resize(1);
        answer.at(1) =  status->giveAlreadyCrack();
        return 1;
    } else if ( type == IST_CrackDirs ) {
        FloatMatrix help;
        status->giveCrackDirs(help);
        answer.resize(9);
        for ( int i = 1; i <= 3; i++ ) {
            answer.at(i) = help.at(1, i);
            answer.at(i + 3) = help.at(2, i);
            answer.at(i + 6) = help.at(3, i);
        }

        return 1;
    } else if ( type == IST_CrackStatuses ) {
        IntArray crackStatus;
        answer.resize(3);
        status->giveCrackStatus(crackStatus);
        for ( int i = 1; i <= 3; i++ ) {
            answer.at(i) = crackStatus.at(i);
        }

        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}


InternalStateValueType
RCM2Material :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_CrackedFlag ) {
        return ISVT_SCALAR;
    } else if ( type == IST_CrackDirs ) {
        return ISVT_VECTOR;
    } else if ( type == IST_CrackStatuses ) {
        return ISVT_VECTOR;
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }
}


int
RCM2Material :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( type == IST_CrackedFlag ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( type ==  IST_CrackDirs ) {
        answer.resize(9);
        for ( int i = 1; i <= 9; i++ ) {
            answer.at(i) = i;
        }

        return 1;
    } else if ( type == IST_CrackStatuses ) {
        answer.resize(3);
        for ( int i = 1; i <= 3; i++ ) {
            answer.at(i) = i;
        }

        return 1;
    } else {
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}


int
RCM2Material :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( type == IST_CrackedFlag ) {
        return 1;
    } else if ( type == IST_CrackDirs ) {
        return 9;
    } else if ( type == IST_CrackStatuses ) {
        return 3;
    } else {
        return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}


void
RCM2Material :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                              MatResponseForm form,
                                              MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *atTime)
{
    //
    // returns receiver 3d material matrix
    //
    this->giveMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
}



void
RCM2Material :: givePlaneStressStiffMtrx(FloatMatrix &answer,
                                         MatResponseForm form,
                                         MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *atTime)

//
// returns receiver's 2dPlaneStressMtrx
// (2dPlaneStres ==> sigma_z = tau_xz = tau_yz = 0.)
//
// standard method from Material Class overloaded, because no inversion is needed.
// the reduction from 3d case will not work
// this implementation should be faster.
{
    this->giveMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
}


void
RCM2Material :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                         MatResponseForm form,
                                         MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *atTime)

//
// return receiver's 2dPlaneStrainMtrx constructed from
// general 3dMatrialStiffnessMatrix
// (2dPlaneStrain ==> eps_z = gamma_xz = gamma_yz = 0.)
//
{
    this->giveMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
}


void
RCM2Material :: give1dStressStiffMtrx(FloatMatrix &answer,
                                      MatResponseForm form,
                                      MatResponseMode mode,
                                      GaussPoint *gp,
                                      TimeStep *atTime)

//
// returns receiver's 1dMaterialStiffnessMAtrix
// (1d case ==> sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.)
{
    this->giveMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
}



void
RCM2Material :: give2dBeamLayerStiffMtrx(FloatMatrix &answer,
                                         MatResponseForm form,
                                         MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *atTime)
//
// returns receiver's 2dBeamLayerStiffMtrx.
// (2dPlaneStres ==> sigma_z = tau_xz = tau_yz = 0.)
//
// standard method from Material Class overloaded, because no inversion is needed.
// the reduction from 3d case will not work
// this implementation should be faster.
{
    this->giveMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
}



void
RCM2Material :: give2dPlateLayerStiffMtrx(FloatMatrix &answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *atTime)
//
// returns receiver's 2dPlateLayerMtrx
// (2dPlaneStres ==> sigma_z = tau_xz = tau_yz = 0.)
//
// standard method from Material Class overloaded, because no inversion is needed.
// the reduction from 3d case will not work
// this implementation should be faster.
{
    this->giveMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
}



void
RCM2Material :: give3dShellLayerStiffMtrx(FloatMatrix &answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *atTime)
//
// returns receiver's 2dPlaneStressMtrx constructed from
// general 3dMatrialStiffnessMatrix
// (2dPlaneStres ==> sigma_z = tau_xz = tau_yz = 0.)
//
// standard method from Material Class overloaded, because no inversion is needed.
// the reduction from 3d case will not work
// this implementation should be faster.
{
    this->give2dPlateLayerStiffMtrx(answer, form, mode, gp, atTime);
}








RCM2MaterialStatus :: RCM2MaterialStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), crackStatuses(3), tempCrackStatuses(3),
    maxCrackStrains(3), tempMaxCrackStrains(3), crackStrainVector(3),
    oldCrackStrainVector(3), crackDirs(3, 3), tempCrackDirs(3, 3), charLengths(3),
    //minEffStrainsForFullyOpenCrack(3),
    principalStrain(3), oldPrincipalStrain(3),
    principalStress(3), oldPrincipalStress(3), crackMap(3) {
    for ( int i = 1; i <= 3; i++ ) {
        crackDirs.at(i, i) = tempCrackDirs.at(i, i) = 1.0;
    }
}


RCM2MaterialStatus :: ~RCM2MaterialStatus()
{ }



int
RCM2MaterialStatus :: isCrackActive(int i) const
//
// return nonzero if crack is active
//
{
    if ( crackStrainVector.at(i) > 0. ) {
        return 1;
    }

    // handle situation for new crack with crackStrainVector->at(i) = 0.
    if ( crackMap.at(i) ) {
        return 1;
    }

    return 0;
}

void
RCM2MaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    int i;
    char s [ 11 ];

    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
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

    fprintf(file, "}\n");
}

int
RCM2MaterialStatus :: giveNumberOfActiveCracks() const
//
// return number of active cracks, computed upon CrackStatus
//
{
    int i, answer = 0;

    for ( i = 1; i <= 3; i++ ) {
        if ( oldCrackStrainVector.at(i) > 0. ) {
            answer++;
        }
    }

    return answer;
}


int
RCM2MaterialStatus :: giveNumberOfTempActiveCracks() const
//
// return number of active cracks, computed upon tempCrackStatus
//
{
    int i, answer = 0;

    for ( i = 1; i <= 3; i++ ) {
        if ( isCrackActive(i) ) {
            answer++;
        }
    }

    return answer;
}



void
RCM2MaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    int i;

    StructuralMaterialStatus :: initTempStatus();

    tempCrackStatuses = crackStatuses;
    tempMaxCrackStrains = maxCrackStrains;
    tempCrackDirs = crackDirs;
    crackStrainVector = oldCrackStrainVector;

    principalStrain = oldPrincipalStrain;
    principalStress = oldPrincipalStress;

    for ( i = 1; i <= 3; i++ ) {
        crackMap.at(i) = 0;
    }
}



void
RCM2MaterialStatus :: updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    StructuralMaterialStatus :: updateYourself(atTime);

    crackStatuses = tempCrackStatuses;
    maxCrackStrains = tempMaxCrackStrains;
    crackDirs = tempCrackDirs;
    oldCrackStrainVector = crackStrainVector;

    oldPrincipalStrain = principalStrain;
    oldPrincipalStress = principalStress;
}



contextIOResultType
RCM2MaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data

    if ( ( iores = crackStatuses.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = maxCrackStrains.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackDirs.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    //if ((iores = minEffStrainsForFullyOpenCrack.storeYourself(stream)) != CIO_OK) THROW_CIOERR(iores);
    if ( ( iores = charLengths.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackStrainVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = oldCrackStrainVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = oldPrincipalStrain.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = principalStrain.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = oldPrincipalStress.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = principalStress.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackMap.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

contextIOResultType
RCM2MaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data

    if ( ( iores = crackStatuses.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = maxCrackStrains.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackDirs.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    //if ((iores = minEffStrainsForFullyOpenCrack.restoreYourself(stream)) != CIO_OK) THROW_CIOERR(iores);
    if ( ( iores = charLengths.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackStrainVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = oldCrackStrainVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = oldPrincipalStrain.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = principalStrain.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = oldPrincipalStress.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = principalStress.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackMap.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK; // return succes
}
} // end namespace oofem
