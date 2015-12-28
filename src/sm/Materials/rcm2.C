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

#include "rcm2.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "contextioerr.h"
#include "mathfem.h"

#include <cstring>

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
    return mode == _3dMat || mode == _PlaneStress ||
           mode == _PlaneStrain || mode == _1dMat ||
           mode == _PlateLayer || mode == _2dBeamLayer;
}


/*
 * void
 * RCM2Material :: computeTrialStressIncrement (FloatArray& answer, GaussPoint *gp,
 *                     const FloatArray &strainIncrement,
 *                          TimeStep* tStep)
 * //
 * // returns trial stress according to strainIncrement
 * //
 * {
 * FloatMatrix def;
 *
 * this -> giveEffectiveMaterialStiffnessMatrix (def, TangentStiffness,gp,tStep);
 * answer.beProductOf (def, strainIncrement);
 * }
 */


void
RCM2Material :: giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
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
    this->giveEffectiveMaterialStiffnessMatrix(answer, mode, gp, tStep);
    // def is full matrix for current gp stress strain mode.
}


void
RCM2Material :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                     const FloatArray &totalStrain,
                                     TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    FloatArray princStress, reducedStressVector, reducedAnswer;
    FloatArray reducedTotalStrainVector, principalStrain, strainVector;
    FloatMatrix tempCrackDirs;
    RCM2MaterialStatus *status = static_cast< RCM2MaterialStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, tStep, VM_Total);
    //

    StructuralMaterial :: giveFullSymVectorForm( strainVector, reducedTotalStrainVector, gp->giveMaterialMode() );

    tempCrackDirs = status->giveTempCrackDirs();
    this->computePrincipalValDir(principalStrain, tempCrackDirs,
                                 strainVector,
                                 principal_strain);

    this->giveRealPrincipalStressVector3d(princStress, gp, principalStrain, tempCrackDirs, tStep);
    princStress.resizeWithValues(6);

    tempCrackDirs = status->giveTempCrackDirs();
    this->transformStressVectorTo(answer, tempCrackDirs, princStress, 1);

    status->letTempStrainVectorBe(totalStrain);
    StructuralMaterial :: giveReducedSymVectorForm( reducedStressVector, answer, gp->giveMaterialMode() );
    status->letTempStressVectorBe(reducedStressVector);
    this->updateCrackStatus(gp, status->giveCrackStrainVector());

    StructuralMaterial :: giveReducedSymVectorForm( reducedAnswer, answer, gp->giveMaterialMode() );
    answer = reducedAnswer;
}


void
RCM2Material :: giveRealPrincipalStressVector3d(FloatArray &answer, GaussPoint *gp,
                                                FloatArray &principalStrain,
                                                FloatMatrix &tempCrackDirs,
                                                TimeStep *tStep)
//
// returns real principal stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
// updates principal strain and stress of the receiver's status.
//
{
    int ind;
    double maxErr;
    FloatArray crackStrainVector;
    FloatArray strainIncrement, crackStrainIterativeIncrement;
    FloatArray prevPrincipalStrain;
    FloatArray dSigma;
    FloatArray elastStrain, sigmaEl, sigmaCr(3);
    FloatArray fullDSigma;
    IntArray activatedCracks, crackMapping;
    FloatMatrix dcr, de, decr, fullDecr, crackDirs;
    RCM2MaterialStatus *status = static_cast< RCM2MaterialStatus * >( this->giveStatus(gp) );

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
    crackStrainVector = status->giveCrackStrainVector(); // local one
    crackDirs = status->giveCrackDirs();
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
    prevPrincipalStrain = status->givePrevPrincStrainVector();
    strainIncrement.beDifferenceOf(principalStrain, prevPrincipalStrain);
    status->letPrincipalStrainVectorBe(principalStrain);

    this->giveNormalElasticStiffnessMatrix(de, false, TangentStiffness,
                                           gp, tStep, tempCrackDirs);
    //
    // construct mapping matrix of active cracks
    // this mapping will dynamically change as
    // some crack can unlo or reload
    //
    this->updateActiveCrackMap(gp);
    crackMapping = status->giveCrackMap();
    // start iteration until stress computed from elastic increment
    // is equal to stress computed from cracking strain increment
    // we do this computation in reduced stress strain space
    dSigma.clear();
    for ( int iter = 1; iter <= 20; iter++ ) {
        //
        // first check if already cracked
        //
        if ( status->giveNumberOfTempActiveCracks() ) {
            // active crack exist
            this->giveCrackedStiffnessMatrix(dcr, TangentStiffness, gp, tStep);
            fullDecr = de;
            fullDecr.add(dcr);
            decr.resize( crackMapping.maximum(), crackMapping.maximum() );
            decr.zero();
            decr.assemble(fullDecr, crackMapping, crackMapping);

            if ( dSigma.giveSize() == 0 ) {
                fullDSigma.beProductOf(de, strainIncrement);
                dSigma.resize( crackMapping.maximum() );
                dSigma.zero();
                dSigma.assemble(fullDSigma, crackMapping);
                //dSigma.beSubArrayOf(fullDSigma, crackMapping); //@todo fix! used old version of method beSubArrayOf
            }

            decr.solveForRhs(dSigma, crackStrainIterativeIncrement);
            for ( int i = 1; i <= 3; i++ ) {
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
            for ( int i = 1; i <= 3; i++ ) {
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
            crackMapping = status->giveCrackMap(); // update crackMap
        }

        //

        // compute unbalanced stress
        // dSigma = sigmaEl - sigmaCr for active cracks
        fullDSigma = sigmaEl;
        fullDSigma.subtract(sigmaCr);
        //fullDSigma.printYourself();
        //crackMapping.printYourself();

        dSigma.resize( crackMapping.maximum() );
        dSigma.zero();
        dSigma.assemble(fullDSigma, crackMapping);
        //dSigma.beSubArrayOf(fullDSigma, crackMapping); //@todo fix! used old version of method beSubArrayOf
        if ( dSigma.giveSize() ) {
            // fullDSigma.printYourself();
            //crackMapping.printYourself();
            //dSigma.printYourself();
        }
        //dSigma.printYourself();

        // find max error in dSigma
        // if max err < allovedErr -> stop iteration
        // allowed Err is computed relative to Ft;

        // check only for active cracks
        maxErr = 0.;
        for ( int i = 1; i <= dSigma.giveSize(); i++ ) {
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
    OOFEM_ERROR("convergence not reached");
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
    int upd, activationFlag = 0;
    RCM2MaterialStatus *status = static_cast< RCM2MaterialStatus * >( this->giveStatus(gp) );
    FloatArray localStress;
    FloatMatrix tempCrackDirs;
    const IntArray &crackMap = status->giveCrackMap();

    answer.resize(3);
    answer.zero();
    localStress = princStressVector;
    //
    // local stress is updated according to reached local crack strain
    //
    for ( int i = 1; i <= 3; i++ ) { // loop over each possible crack plane
        // test previous status of each possible crack plane
        upd = 0;
        IntArray indx;
        StructuralMaterial :: giveVoigtSymVectorMask( indx, gp->giveMaterialMode() );
        if ( crackMap.at(i) == 0 && indx.contains(i) ) {
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
                tempCrackDirs = status->giveTempCrackDirs();
                for ( int j = 1; j <= 3; j++ ) {
                    crackPlaneNormal.at(j) = tempCrackDirs.at(j, i);
                }

                Le = this->giveCharacteristicElementLength(gp, crackPlaneNormal);
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

    answer.clear();
}


void
RCM2Material :: updateStatusForNewCrack(GaussPoint *gp, int i, double Le)
//
// updates gp status when new crack-plane i is formed with charLength Le
// updates Le and computes and sets minEffStrainForFullyOpenCrack
//
{
    RCM2MaterialStatus *status = static_cast< RCM2MaterialStatus * >( this->giveStatus(gp) );

    if ( Le <= 0 ) {
        OOFEM_ERROR("Element %d returned zero char length", gp->giveElement()->giveNumber() );
    }

    status->setCharLength(i, Le);
    //status ->  setMinCrackStrainsForFullyOpenCrack (i, this->giveMinCrackStrainsForFullyOpenCrack(gp,i));
}

double 
RCM2Material :: giveCharacteristicElementLength(GaussPoint *gp, const FloatArray &crackPlaneNormal) 
{ 
  return gp->giveElement()->giveCharacteristicLength(crackPlaneNormal); 
}

void
RCM2Material :: updateCrackStatus(GaussPoint *gp, const FloatArray &crackStrain)
//
// updates gp records and MatStatus due to cracking.
//
// Updates MatStatus (its temporary variables) to current
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
    double minCrackStrainsForFullyOpenCrack;
    RCM2MaterialStatus *status = static_cast< RCM2MaterialStatus * >( this->giveStatus(gp) );
    const IntArray &crackMap = status->giveCrackMap();

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
    IntArray indx;
    StructuralMaterial :: giveVoigtSymVectorMask( indx, gp->giveMaterialMode() );
    for ( int i = 1; i <= 3; i++ ) { // loop over each possible crack plane
        if ( crackMap.at(i) != 0 && indx.contains(i) ) {
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
    RCM2MaterialStatus *status = static_cast< RCM2MaterialStatus * >( this->giveStatus(gp) );
    int isClosing = 0;

    for ( int i = 1; i <= 3; i++ ) {
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

    crackMap = status->giveCrackMap(); // update crack Map
}


void
RCM2Material :: giveNormalElasticStiffnessMatrix(FloatMatrix &answer,
                                                 bool reduce, MatResponseMode rMode,
                                                 GaussPoint *gp, TimeStep *tStep,
                                                 const FloatMatrix &dir)
//
// return Elastic Stiffness matrix for normal Stresses
// taking into account the directions of normal stresses
// (not supported now)
//
{
    StructuralMaterial *lMat = this->giveLinearElasticMaterial();
    FloatMatrix de, fullAnswer(3, 3);
    IntArray mask;
    int sd;

    FloatMatrix stiff;
    lMat->giveStiffnessMatrix(stiff, rMode, gp, tStep);
    this->giveFullSymMatrixForm( de, stiff, gp->giveMaterialMode() );

    // copy first 3x3 submatrix to answer
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            fullAnswer.at(i, j) = de.at(i, j);
        }
    }

    if ( !reduce ) { // 3x3 full form required
        answer = fullAnswer;
    } else {
        // reduced form for only
        // nonzero normal stresses
        //
        // first find spatial dimension of problem
        sd = 0;

        IntArray indx;
        StructuralMaterial :: giveVoigtSymVectorMask( indx, gp->giveMaterialMode() );
        for ( int i = 1; i <= 3; i++ ) {
            if ( indx.contains(i) ) {
                sd++;
            }
        }

        answer.resize(sd, sd);
        answer.zero();

        // copy fullAnswer to reduced one
        StructuralMaterial :: giveInvertedVoigtVectorMask( mask, gp->giveMaterialMode() );
        for ( int i = 1; i <= sd; i++ ) {
            int iidx = mask.findFirstIndexOf(i);
            for ( int j = 1; j <= sd; j++ ) {
                int jidx = mask.findFirstIndexOf(j);
                if ( iidx && jidx ) {
                    answer.at(i, j) = fullAnswer.at(iidx, jidx);
                }
            }
        }
    }
}


void
RCM2Material :: giveEffectiveMaterialStiffnessMatrix(FloatMatrix &answer,
                                                     MatResponseMode rMode, GaussPoint *gp,
                                                     TimeStep *tStep)
//
// returns effective material stiffness matrix in full form
// for gp stress strain mode
//
{
    RCM2MaterialStatus *status = static_cast< RCM2MaterialStatus * >( this->giveStatus(gp) );
    StructuralMaterial *lMat = static_cast< StructuralMaterial * >( this->giveLinearElasticMaterial() );
    int numberOfActiveCracks = status->giveNumberOfTempActiveCracks();
    int indi, indj, ii, jj;
    double G, princStressDis, princStrainDis;
    FloatMatrix de, invDe, compliance, dcr, d, df, t, tempCrackDirs;
    FloatArray principalStressVector, principalStrainVector;
    IntArray mask;

    if ( ( rMode == ElasticStiffness ) || ( numberOfActiveCracks == 0 ) ) {
        lMat->giveStiffnessMatrix(answer, rMode, gp, tStep);
        return;
    }

    // this->updateActiveCrackMap(gp) must be done after restart.
    this->updateActiveCrackMap(gp);
    tempCrackDirs = status->giveTempCrackDirs();
    this->giveNormalElasticStiffnessMatrix(de, true, rMode, gp, tStep,
                                           tempCrackDirs);
    invDe.beInverseOf(de);
    this->giveCrackedStiffnessMatrix(dcr, rMode, gp, tStep);
    StructuralMaterial :: giveVoigtSymVectorMask( mask, gp->giveMaterialMode() );
    compliance.resize( mask.giveSize(), mask.giveSize() );

    // we will set
    // first we set compliances for normal streses in
    // local coordinate system defined by crackplane
    for ( int i = 1; i <= 3; i++ ) {
        IntArray indx;
        StructuralMaterial :: giveVoigtSymVectorMask( indx, gp->giveMaterialMode() );
        if ( ( indi = indx.findFirstIndexOf(i) ) ) {
            for ( int j = 1; j <= 3; j++ ) {
                if ( ( indj = indx.findFirstIndexOf(j) ) ) {
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

    principalStressVector = status->getPrincipalStressVector();
    principalStrainVector = status->getPrincipalStrainVector();

    // now remain to set shears
    G = this->give(pscm_G, gp);

    IntArray indx;
    StructuralMaterial :: giveVoigtSymVectorMask( indx, gp->giveMaterialMode() );
    for ( int i = 4; i <= 6; i++ ) {
        if ( ( indi = indx.findFirstIndexOf(i) ) ) {
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
    StructuralMaterial :: giveVoigtSymVectorMask( mask, gp->giveMaterialMode() );
    df.resize(6, 6);
    df.zero();
    df.assemble(d, mask);
    //
    // final step - transform stiffnes to global c.s
    //

    this->giveStressVectorTranformationMtrx(t, tempCrackDirs, 1);
    df.rotatedWith(t, 't');

    StructuralMaterial :: giveReducedSymMatrixForm( answer, df, gp->giveMaterialMode() );
}


void
RCM2Material :: giveCrackedStiffnessMatrix(FloatMatrix &answer,
                                           MatResponseMode rMode,
                                           GaussPoint *gp,
                                           TimeStep *tStep)
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
    RCM2MaterialStatus *status = static_cast< RCM2MaterialStatus * >( this->giveStatus(gp) );
    int numberOfActiveCracks = status->giveNumberOfTempActiveCracks();
    const IntArray &crackMap = status->giveCrackMap();
    if ( numberOfActiveCracks == 0 ) {
        answer.clear();
        return;
    }

    answer.resize(3, 3);
    answer.zero();

    // loop over each active crack plane
    for ( int i = 1; i <= 3; i++ ) {
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
    int indx = 1;
    RCM2MaterialStatus *status = static_cast< RCM2MaterialStatus * >( this->giveStatus(gp) );
    IntArray crackMap = status->giveCrackMap();

    //if (crackMap == NULL) OOFEM_ERROR("NULL pointer encountered");

    for ( int i = 1; i <= 3; i++ ) {
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
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, Gf, _IFT_RCM2Material_gf);
    IR_GIVE_FIELD(ir, Ft, _IFT_RCM2Material_ft);

    return this->giveLinearElasticMaterial()->initializeFrom(ir);
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

    if ( propertyDictionary.includes(aProperty) ) {
        value = propertyDictionary.at(aProperty);
    } else {
        if ( linearElasticMaterial ) {
            value = this->linearElasticMaterial->give(aProperty, gp);
        } else {
            OOFEM_ERROR("property not defined");
        }
    }

    return value;
}


int
RCM2Material :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    RCM2MaterialStatus *status = static_cast< RCM2MaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_CrackedFlag ) {
        answer.resize(1);
        answer.at(1) =  status->giveAlreadyCrack();
        return 1;
    } else if ( type == IST_CrackDirs ) {
        const FloatMatrix &help = status->giveCrackDirs();
        answer.resize(9);
        for ( int i = 1; i <= 3; i++ ) {
            answer.at(i) = help.at(1, i);
            answer.at(i + 3) = help.at(2, i);
            answer.at(i + 6) = help.at(3, i);
        }

        return 1;
    } else if ( type == IST_CrackStatuses ) {
        const IntArray &crackStatus = status->giveCrackStatus();
        answer.resize(3);
        for ( int i = 1; i <= 3; i++ ) {
            answer.at(i) = crackStatus.at(i);
        }

        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

void
RCM2Material :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                              MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *tStep)
{
    //
    // returns receiver 3d material matrix
    //
    this->giveMaterialStiffnessMatrix(answer, mode, gp, tStep);
}


void
RCM2Material :: givePlaneStressStiffMtrx(FloatMatrix &answer,
                                         MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *tStep)

//
// returns receiver's 2dPlaneStressMtrx
// (2dPlaneStres ==> sigma_z = tau_xz = tau_yz = 0.)
//
// standard method from Material Class overloaded, because no inversion is needed.
// the reduction from 3d case will not work
// this implementation should be faster.
{
    this->giveMaterialStiffnessMatrix(answer, mode, gp, tStep);
}


void
RCM2Material :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                         MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *tStep)

//
// return receiver's 2dPlaneStrainMtrx constructed from
// general 3dMatrialStiffnessMatrix
// (2dPlaneStrain ==> eps_z = gamma_xz = gamma_yz = 0.)
//
{
    this->giveMaterialStiffnessMatrix(answer, mode, gp, tStep);
}


void
RCM2Material :: give1dStressStiffMtrx(FloatMatrix &answer,
                                      MatResponseMode mode,
                                      GaussPoint *gp,
                                      TimeStep *tStep)

//
// returns receiver's 1dMaterialStiffnessMAtrix
// (1d case ==> sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.)
{
    this->giveMaterialStiffnessMatrix(answer, mode, gp, tStep);
}


void
RCM2Material :: give2dBeamLayerStiffMtrx(FloatMatrix &answer,
                                         MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *tStep)
//
// returns receiver's 2dBeamLayerStiffMtrx.
// (2dPlaneStres ==> sigma_z = tau_xz = tau_yz = 0.)
//
// standard method from Material Class overloaded, because no inversion is needed.
// the reduction from 3d case will not work
// this implementation should be faster.
{
    this->giveMaterialStiffnessMatrix(answer, mode, gp, tStep);
}


void
RCM2Material :: givePlateLayerStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode mode,
                                        GaussPoint *gp,
                                        TimeStep *tStep)
//
// returns receiver's 2dPlateLayerMtrx
// (2dPlaneStres ==> sigma_z = tau_xz = tau_yz = 0.)
//
// standard method from Material Class overloaded, because no inversion is needed.
// the reduction from 3d case will not work
// this implementation should be faster.
{
    this->giveMaterialStiffnessMatrix(answer, mode, gp, tStep);
}




RCM2MaterialStatus :: RCM2MaterialStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), crackStatuses(3), tempCrackStatuses(3),
    maxCrackStrains(3), tempMaxCrackStrains(3), crackStrainVector(3),
    oldCrackStrainVector(3), crackDirs(3, 3), tempCrackDirs(3, 3), charLengths(3),
    //minEffStrainsForFullyOpenCrack(3),
    principalStrain(3), oldPrincipalStrain(3),
    principalStress(3), oldPrincipalStress(3), crackMap(3)
{
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
    int answer = 0;

    for ( int i = 1; i <= 3; i++ ) {
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
    int answer = 0;

    for ( int i = 1; i <= 3; i++ ) {
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
    StructuralMaterialStatus :: initTempStatus();

    tempCrackStatuses = crackStatuses;
    tempMaxCrackStrains = maxCrackStrains;
    tempCrackDirs = crackDirs;
    crackStrainVector = oldCrackStrainVector;

    principalStrain = oldPrincipalStrain;
    principalStress = oldPrincipalStress;

    for ( int i = 1; i <= 3; i++ ) {
        crackMap.at(i) = 0;
    }
}



void
RCM2MaterialStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    StructuralMaterialStatus :: updateYourself(tStep);

    crackStatuses = tempCrackStatuses;
    maxCrackStrains = tempMaxCrackStrains;
    crackDirs = tempCrackDirs;
    oldCrackStrainVector = crackStrainVector;

    oldPrincipalStrain = principalStrain;
    oldPrincipalStress = principalStress;
}



contextIOResultType
RCM2MaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
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

    if ( ( iores = crackStatuses.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = maxCrackStrains.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackDirs.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    //if ((iores = minEffStrainsForFullyOpenCrack.storeYourself(stream)) != CIO_OK) THROW_CIOERR(iores);
    if ( ( iores = charLengths.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackStrainVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = oldCrackStrainVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = oldPrincipalStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = principalStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = oldPrincipalStress.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = principalStress.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackMap.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

contextIOResultType
RCM2MaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
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

    if ( ( iores = crackStatuses.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = maxCrackStrains.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackDirs.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    //if ((iores = minEffStrainsForFullyOpenCrack.restoreYourself(stream)) != CIO_OK) THROW_CIOERR(iores);
    if ( ( iores = charLengths.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackStrainVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = oldCrackStrainVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = oldPrincipalStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = principalStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = oldPrincipalStress.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = principalStress.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackMap.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK; // return succes
}
} // end namespace oofem
