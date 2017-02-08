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

#include "fcm.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "stressvector.h"

#include <cstring>

namespace oofem {
FCMMaterial :: FCMMaterial(int n, Domain *d) : StructuralMaterial(n, d)
    //
    // constructor
    //
{
    ecsMethod = ECSM_Unknown;
    linearElasticMaterial = NULL;
}


FCMMaterial :: ~FCMMaterial()
//
// destructor
//
{}

int
FCMMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    return mode == _3dMat || mode == _PlaneStress || mode == _PlaneStrain;
}

void
FCMMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                    const FloatArray &totalStrain,
                                    TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    // g = global direction
    // l = local (crack) direction

    // stress and strain transformation matrices
    FloatMatrix epsG2L, sigmaL2G;
    // stiffness matrices (total, elastic, cracking)
    FloatMatrix D, De, Dcr;
    FloatMatrix principalDirs;

    // STRESSES
    FloatArray stressIncrement_l, stressVector_g, trialStress_g, trialStress_l, sigmaResid, sigmaElast_l, sigmaCrack_l;

    // STRAINS
    FloatArray reducedStrain_g, reducedStrain_l, strainIncrement_g, strainIncrement_l, crackStrain, crackStrainIncrement, elasticStrain_l;

    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );
    MaterialMode mMode = gp->giveMaterialMode();

    int nCr = status->giveNumberOfCracks();
    int nMaxCr = status->giveMaxNumberOfCracks(gp);

    double maxErr;
    int iterLimit;
    double maxTau;

    // for the bisection method:
    int index, indexCount;
    bool exitWhileCond = false;
    bool cancelledStrainFlag = false;

    bool plus, minus;
    FloatArray crackStrainPlus, crackStrainMinus, residPlus, residMinus; //local

    index = indexCount = 0;

    this->initTempStatus(gp);

    // get elastic stiffness matrix
    this->giveStiffnessMatrix(De, ElasticStiffness, gp, tStep);


    // NOW INDEPENDENTLY OF THE HISTORY GET PRINCIPAL STRESSES / STRESSES IN THE CRACK DIRS

    // ELASTIC CASE - NO CRACKING SO FAR (in prevs steps)
    if ( nCr == 0 ) {
        this->giveStressDependentPartOfStrainVector(reducedStrain_g, gp, totalStrain, tStep, VM_Total);
        trialStress_g.beProductOf(De, reducedStrain_g);

        this->computePrincipalValDir(trialStress_l, principalDirs, trialStress_g, principal_stress);

        // REMAINS ELASTIC or element cracking is prevented from above
        if ( !this->isStrengthExceeded( principalDirs, gp, tStep, 1, trialStress_l.at(1) ) ) {
            answer = trialStress_g;
            status->letTempStrainVectorBe(totalStrain);
            status->letTempStressVectorBe(trialStress_g);
            return;

            // STARTS CRACKING - 1st crack
        } else {
            this->initializeCrack(gp, principalDirs, 1);

            for ( int iCrack = 2; iCrack <= min(nMaxCr, this->nAllowedCracks); iCrack++ ) {
                if ( this->isStrengthExceeded( principalDirs, gp, tStep, iCrack, trialStress_l.at(iCrack) ) ) { // for "nonlocal model"
                    this->initializeCrack(gp, principalDirs, iCrack);
                }
            }
        }

        // CRACKING OCCURED IN PREVIOUS STEPS
    } else {
        // check other principal directions - if strength is not exceeded, take transformation matrices from previous analysi
        if ( nCr < nMaxCr ) {
            // equilibrated global stress
            stressVector_g = status->giveStressVector();
            // get strain increment
            this->giveStressDependentPartOfStrainVector(reducedStrain_g, gp, totalStrain, tStep, VM_Incremental);
            strainIncrement_g.beDifferenceOf( reducedStrain_g, status->giveStrainVector() );

            // stress increment in crack direction
            trialStress_g.beProductOf(De, strainIncrement_g);
            trialStress_g.add(stressVector_g);

            principalDirs = status->giveCrackDirs();

            for ( int iCrack = nCr + 1; iCrack <= min(nMaxCr, this->nAllowedCracks); iCrack++ ) {
                // condition to guarantee that third  crack won't be initiated if second crack does not exist
                if ( ( status->giveNumberOfTempCracks() + 1 ) < iCrack ) {
                    break;
                }

                if ( !this->checkStrengthCriterion(principalDirs, trialStress_g, gp, tStep, iCrack) ) { // = strength of composite
                    // second and third crack is now initialized
                    this->initializeCrack(gp, principalDirs, iCrack);
                }
            } // end for
        } // end nMaxCr
    } // cracking in previous steps

    sigmaL2G = status->giveL2GStressVectorTransformationMtrx();
    epsG2L = status->giveG2LStrainVectorTransformationMtrx();


    // PREPARE STRAIN AND STRESS INCREMENT IN CRACK DIRECTIONS

    if ( ( nCr == 0 ) || ( nCr == nMaxCr ) ) { // not to do the same job twice
        // get strain increment
        this->giveStressDependentPartOfStrainVector(reducedStrain_g, gp, totalStrain, tStep, VM_Incremental);
        strainIncrement_g.beDifferenceOf( reducedStrain_g, status->giveStrainVector() );
    }


    // strain and strain increment in crack direction
    strainIncrement_l = strainIncrement_g;
    strainIncrement_l.rotatedWith(epsG2L, 'n'); // from global to local

    reducedStrain_l = reducedStrain_g;
    reducedStrain_l.rotatedWith(epsG2L, 'n'); // from global to local

    // stress increment in crack direction
    stressIncrement_l.beProductOf(De, strainIncrement_l);


    // SIMILARLY TO RCM2: THE STRESS INCREMENT IN ELASTIC
    // AND CRACKING UNIT MUST BE EQUAL - ITERATE

    crackStrain = status->giveCrackStrainVector();

    iterLimit = 20;

    for ( int iter = 1; iter <= iterLimit; iter++ ) {
        this->giveLocalCrackedStiffnessMatrix(Dcr, TangentStiffness, gp, tStep);

        if ( iter == 1 ) {
            // residuum
            sigmaResid = stressIncrement_l;
        }

        D = De;
        D.add(Dcr);

        D.solveForRhs(sigmaResid, crackStrainIncrement);

        // update and store cracking strain
        crackStrain.add(crackStrainIncrement);

        status->setTempCrackStrainVector(crackStrain);

        // update statuses
        this->updateCrackStatus(gp);

        // compute and compare stress in elastic and cracking unit:
        // EL
        elasticStrain_l = reducedStrain_l;
        elasticStrain_l.subtract(crackStrain);
        sigmaElast_l.beProductOf(De, elasticStrain_l);
        // CR
        this->giveLocalCrackedStiffnessMatrix(Dcr, SecantStiffness, gp, tStep);
        //   sigmaCrack_l.beProductOf(Dr, crackStrain); // matrix is diagonal
        sigmaCrack_l.resize( crackStrain.giveSize() );
        for ( int i = 1; i <= crackStrain.giveSize(); i++ ) {
            sigmaCrack_l.at(i) = Dcr.at(i, i) * crackStrain.at(i);

            switch ( mMode ) {
            case  _PlaneStress:

                if ( i == 3 ) {
                    maxTau = this->maxShearStress(gp, 6);

                    if ( sigmaCrack_l.at(i) > maxTau ) {
                        sigmaCrack_l.at(i) = maxTau;
                    }
                }
                break;

            case _3dMat:

                if ( i >= 4 ) {
                    maxTau = this->maxShearStress(gp, i);

                    if ( sigmaCrack_l.at(i) > maxTau ) {
                        sigmaCrack_l.at(i) = maxTau;
                    }
                }
                break;

            case _PlaneStrain: // check

                if ( i == 4 ) {
                    maxTau = this->maxShearStress(gp, 6);

                    if ( sigmaCrack_l.at(i) > maxTau ) {
                        sigmaCrack_l.at(i) = maxTau;
                    }
                }
                break;

            default:
                OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
            }
        }

        // residuum
        sigmaResid = sigmaElast_l;
        sigmaResid.subtract(sigmaCrack_l);


        maxErr = 0.;
        for ( int i = 1; i <= sigmaResid.giveSize(); i++ ) {
            if ( fabs( sigmaResid.at(i) ) > maxErr ) {
                maxErr = fabs( sigmaResid.at(i) );
            }
        }


        if  ( maxErr < fcm_TOLERANCE * this->giveTensileStrength(gp) ) {
            break;
        }
    }


    // modified iteration method - first gradient and then bisection method

    while ( maxErr > fcm_TOLERANCE * this->giveTensileStrength(gp) ) {
        for ( int i = 1; i <= sigmaResid.giveSize(); i++ ) {
            if ( fabs( sigmaResid.at(i) ) == maxErr ) {
                // found the same index as before -> do not do the same job over, exit while
                if ( index == i ) {
                    exitWhileCond = true;
                    break;
                } else {
                    // the crack strain index is not the same as before - clear the following flag
                    cancelledStrainFlag = false;
                }

                index = i;

                break;
            }
        }

        if ( exitWhileCond ) {
            break;
        }

        iterLimit = 100;

        plus = false;
        minus = false;

        crackStrainPlus.zero();
        crackStrainMinus.zero();
        residPlus.zero();
        residMinus.zero();

        indexCount++;

        if ( indexCount > 10 ) {
            OOFEM_WARNING("Fixed crack model: Local equilibrium not reached!, max. stress error %f", maxErr);
            break;
        }


        for ( int iter = 1; iter <= iterLimit; iter++ ) {
            if ( iter == iterLimit ) {
                OOFEM_WARNING("Fixed crack model: Local equilibrium not reached!, max. stress error %f", maxErr);
            }


            if ( ( iter == 1 ) && ( indexCount == 1 ) ) {
                crackStrain.zero();
                crackStrain.add( status->giveCrackStrainVector() );
                status->setTempCrackStrainVector(crackStrain);
                //sigmaResid.zero();
                sigmaResid = stressIncrement_l;
                this->updateCrackStatus(gp);
            }

            // DO THE FOLLOWING SECTION ONLY UNTIL "PLUS" AND "MINUS" CRACK STRAIN IS FOUND

            if ( ( !plus ) || ( !minus ) ) {
                this->giveLocalCrackedStiffnessMatrix(Dcr, TangentStiffness, gp, tStep);

                D = De;
                D.add(Dcr);

                D.solveForRhs(sigmaResid, crackStrainIncrement);

                // update and store cracking strain
                crackStrain.add(crackStrainIncrement);

                status->setTempCrackStrainVector(crackStrain);

                // update statuses
                this->updateCrackStatus(gp);

                // compute and compare stress in elastic and cracking unit:
                // EL
                elasticStrain_l = reducedStrain_l;
                elasticStrain_l.subtract(crackStrain);
                sigmaElast_l.beProductOf(De, elasticStrain_l);
                // CR
                this->giveLocalCrackedStiffnessMatrix(Dcr, SecantStiffness, gp, tStep);
                //   sigmaCrack_l.beProductOf(Dr, crackStrain); // matrix is diagonal
                sigmaCrack_l.resize( crackStrain.giveSize() );
                for ( int i = 1; i <= crackStrain.giveSize(); i++ ) {
                    sigmaCrack_l.at(i) = Dcr.at(i, i) * crackStrain.at(i);

                    switch ( mMode ) {
                    case  _PlaneStress:

                        if ( i == 3 ) {
                            maxTau = this->maxShearStress(gp, 6);

                            if ( sigmaCrack_l.at(i) > maxTau ) {
                                sigmaCrack_l.at(i) = maxTau;
                            }
                        }
                        break;

                    case _3dMat:

                        if ( i >= 4 ) {
                            maxTau = this->maxShearStress(gp, i);

                            if ( sigmaCrack_l.at(i) > maxTau ) {
                                sigmaCrack_l.at(i) = maxTau;
                            }
                        }
                        break;

                    case _PlaneStrain: // check

                        if ( i == 4 ) {
                            maxTau = this->maxShearStress(gp, 6);

                            if ( sigmaCrack_l.at(i) > maxTau ) {
                                sigmaCrack_l.at(i) = maxTau;
                            }
                        }
                        break;

                    default:
                        OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
                    }
                }

                // residuum
                sigmaResid = sigmaElast_l;
                sigmaResid.subtract(sigmaCrack_l);

                if ( sigmaResid.at(index) > 0 ) {
                    if ( !plus ) { // first positive sigma resid
                        residPlus = sigmaResid;
                    } else {
                        if ( sigmaResid.at(index) < residPlus.at(index) ) { // smaller error
                            residPlus = sigmaResid;
                        }
                    }

                    crackStrainPlus = crackStrain;
                    plus = true;
                } else { // sigmaResid < 0
                    if ( !minus ) { // first negative sigma resid
                        residMinus = sigmaResid;
                    } else {
                        if ( sigmaResid.at(index) > residMinus.at(index) ) { // smaller error
                            residMinus = sigmaResid;
                        }
                    }

                    crackStrainMinus = crackStrain;
                    minus = true;
                }


                maxErr = 0.;
                for ( int i = 1; i <= sigmaResid.giveSize(); i++ ) {
                    if ( fabs( sigmaResid.at(i) ) > maxErr ) {
                        maxErr = fabs( sigmaResid.at(i) );
                    }
                }


                if  ( fabs( sigmaResid.at(index) ) < fcm_TOLERANCE * this->giveTensileStrength(gp) ) {
                    break;
                }
            } else { // we have "plus" and "minus" cracking and the bisection method can start
                crackStrain.zero();
                crackStrain.add(crackStrainMinus);
                crackStrain.add(crackStrainPlus);
                crackStrain.times(0.5);

                status->setTempCrackStrainVector(crackStrain);
                // update statuses - to have correct stiffnesses (unlo-relo vs. softening)
                this->updateCrackStatus(gp);



                // compute and compare stress in elastic and cracking unit:
                // EL
                elasticStrain_l = reducedStrain_l;
                elasticStrain_l.subtract(crackStrain);
                sigmaElast_l.beProductOf(De, elasticStrain_l);
                // CR
                this->giveLocalCrackedStiffnessMatrix(Dcr, SecantStiffness, gp, tStep);
                //   sigmaCrack_l.beProductOf(Dr, crackStrain); // matrix is diagonal
                sigmaCrack_l.resize( crackStrain.giveSize() );
                for ( int i = 1; i <= crackStrain.giveSize(); i++ ) {
                    sigmaCrack_l.at(i) = Dcr.at(i, i) * crackStrain.at(i);

                    switch ( mMode ) {
                    case  _PlaneStress:

                        if ( i == 3 ) {
                            maxTau = this->maxShearStress(gp, 6);

                            if ( sigmaCrack_l.at(i) > maxTau ) {
                                sigmaCrack_l.at(i) = maxTau;
                            }
                        }
                        break;

                    case _3dMat:

                        if ( i >= 4 ) {
                            maxTau = this->maxShearStress(gp, i);

                            if ( sigmaCrack_l.at(i) > maxTau ) {
                                sigmaCrack_l.at(i) = maxTau;
                            }
                        }
                        break;

                    case _PlaneStrain: // check

                        if ( i == 4 ) {
                            maxTau = this->maxShearStress(gp, 6);

                            if ( sigmaCrack_l.at(i) > maxTau ) {
                                sigmaCrack_l.at(i) = maxTau;
                            }
                        }
                        break;

                    default:
                        OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
                    }
                }

                // residuum
                sigmaResid = sigmaElast_l;
                sigmaResid.subtract(sigmaCrack_l);


                maxErr = 0.;
                for ( int i = 1; i <= sigmaResid.giveSize(); i++ ) {
                    if ( fabs( sigmaResid.at(i) ) > maxErr ) {
                        maxErr = fabs( sigmaResid.at(i) );
                    }
                }

                // condition for the given stress component
                if  ( fabs( sigmaResid.at(index) ) < fcm_TOLERANCE * this->giveTensileStrength(gp) ) {
                    break;
                }


                if ( sigmaResid.at(index) > 0 ) {
                    if ( sigmaResid.at(index) < residPlus.at(index) ) { // smaller error
                        residPlus = sigmaResid;
                    }

                    crackStrainPlus = crackStrain;
                } else { // sigmaResid < 0
                    if ( sigmaResid.at(index) > residMinus.at(index) ) { // smaller error
                        residMinus = sigmaResid;
                    }

                    crackStrainMinus = crackStrain;
                }


                if ( fabs( crackStrainMinus.at(index) - crackStrainPlus.at(index) ) < 1.e-18 ) {
                    crackStrain.zero();
                    crackStrain.add(crackStrainMinus);
                    crackStrain.add(crackStrainPlus);
                    crackStrain.times(0.5);

                    //	  if ( fabs (crackStrain.at(index) < 1.e-18 ) ) {
                    crackStrain.at(index) = 0.;
                    cancelledStrainFlag = true;
                    //	  }

#if DEBUG
                    if ( cancelledStrainFlag ) {
                        OOFEM_WARNING("Fixed crack model: cracking strain component %d set to zero", index);
                    }
#endif

                    break;
                }
            } // plus-minus condition
        } // iter looop
    } // maxErr exceeded


    if ( ( maxErr >= fcm_TOLERANCE * this->giveTensileStrength(gp) ) && ( !cancelledStrainFlag ) ) {
        OOFEM_WARNING("Fixed crack model: Local equilibrium not reached!, max. stress error %f", maxErr);
    }

    // set the final-equilibrated crack strain and store the correct corresponding statuses
    status->setTempCrackStrainVector(crackStrain);
    // update statuses - to have correct stiffnesses
    this->updateCrackStatus(gp);


    // correct non-zeros, statuses, etc...
    for ( int i = 1; i <= nMaxCr; i++ ) {
#if DEBUG
        // check for NaNs
        if ( crackStrain.at(i) != crackStrain.at(i) ) {
            OOFEM_ERROR( "crack strain is NaN: %f", crackStrain.at(i) );
        }
#endif

        if ( ( status->giveTempCrackStatus(i) == pscm_NONE ) && ( crackStrain.at(i) <= 0. ) ) {
            status->setTempCrackStrain(i, 0.);
        } else if ( ( status->giveTempCrackStatus(i) == pscm_JUST_INIT ) && ( crackStrain.at(i) <= fcm_SMALL_STRAIN ) ) {
            if ( i + 1 > status->giveMaxNumberOfCracks(gp) ) {
                status->setTempCrackStatus(i, pscm_NONE);
                status->setTempCrackStrain(i, 0.);
                if ( i == 1 ) {
                    principalDirs.zero();
                    status->setCrackDirs(principalDirs);
                }

                // can't change crack if the subsequent crack exists
            } else if ( status->giveTempCrackStatus(i + 1) != pscm_NONE ) {} else {
                status->setTempCrackStatus(i, pscm_NONE);
                status->setTempCrackStrain(i, 0.);
                if ( i == 1 ) {
                    principalDirs.zero();
                    status->setCrackDirs(principalDirs);
                }
            }
        } else if ( crackStrain.at(i) < 0. ) {
            status->setTempCrackStrain(i, 0.);
        }
    }

    // at first crack initiation there should not be any other cracking strains
    if ( ( status->giveTempCrackStatus(1) == pscm_JUST_INIT ) && ( status->giveTempCrackStatus(2) == pscm_NONE ) ) {
        status->setTempCrackStrain(2, 0.);
        status->setTempCrackStrain(3, 0.);
    }


    // convert from local to global coordinates
    sigmaElast_l.rotatedWith(sigmaL2G, 'n');

    answer = sigmaElast_l;
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(sigmaElast_l);

    return;
}


void
FCMMaterial :: initializeCrack(GaussPoint *gp, FloatMatrix &base, int nCrack)
{
    MaterialMode mMode = gp->giveMaterialMode();
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    FloatArray crackVector(3);
    crackVector.zero();

    FloatMatrix epsG2L, sigmaG2L, epsL2G, sigmaL2G;

    int nMaxCracks;
    double Le;

    for ( int i = 1; i <= base.giveNumberOfRows(); i++ ) {
        crackVector.at(i) = base.at(i, nCrack);
    }

    Le = this->giveCharacteristicElementLength(gp, crackVector);
    status->setCharLength(nCrack, Le);

    this->checkSnapBack(gp, nCrack);

    status->setTempCrackStatus(nCrack, pscm_JUST_INIT);
    //status->setTempCrackStatus(nCrack, pscm_SOFTENING);

    // change the base crack vector and transformation matrices
    // when the initiated crack is not the last possible one

    // return for the second crack in plane stress - base is set and transformation matrices are given by the first crack
    if ( ( mMode == _PlaneStress ) && ( nCrack == 2 ) ) {
        return;
    }

    nMaxCracks = status->giveMaxNumberOfCracks(gp);

    if ( nCrack <= nMaxCracks ) {
        // make sure that the matrix of base vectors is right-handed
        if ( nMaxCracks == 3 ) {
            base.at(1, 3) = base.at(2, 1) * base.at(3, 2) - base.at(3, 1) * base.at(2, 2);
            base.at(2, 3) = base.at(3, 1) * base.at(1, 2) - base.at(1, 1) * base.at(3, 2);
            base.at(3, 3) = base.at(1, 1) * base.at(2, 2) - base.at(2, 1) * base.at(1, 2);
        }

        //    status->setTempCrackDirs(base);
        status->setCrackDirs(base);

        if ( mMode == _PlaneStress ) {
            this->givePlaneStressVectorTranformationMtrx(sigmaG2L, base, 0);
            this->give2DStrainVectorTranformationMtrx(epsG2L, base, 0);

            this->givePlaneStressVectorTranformationMtrx(sigmaL2G, base, 1);
            this->give2DStrainVectorTranformationMtrx(epsL2G, base, 1);
        } else {
            this->giveStressVectorTranformationMtrx(sigmaG2L, base, 0);
            this->giveStrainVectorTranformationMtrx(epsG2L, base, 0);

            this->giveStressVectorTranformationMtrx(sigmaL2G, base, 1);
            this->giveStrainVectorTranformationMtrx(epsL2G, base, 1);
        }

        status->setG2LStressVectorTransformationMtrx(sigmaG2L);
        status->setG2LStrainVectorTransformationMtrx(epsG2L);

        status->setL2GStressVectorTransformationMtrx(sigmaL2G);
        status->setL2GStrainVectorTransformationMtrx(epsL2G);
    }
}


bool
FCMMaterial :: isIntact(GaussPoint *gp, int icrack) {
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    if ( icrack >= 4 ) {
        OOFEM_ERROR("Unexpected crack number");
    }

    if ( ( status->giveTempCrackStatus(icrack) != pscm_NONE ) && ( status->giveTempCrackStatus(icrack) != pscm_CLOSED ) ) {
        return false;
    } else {
        return true;
    }
}




bool
FCMMaterial :: isIntactForShear(GaussPoint *gp, int i) {
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    int normal_1, normal_2;

    if ( i == 4 ) { // y-z
        normal_1 = 2;
        normal_2 = 3;
    } else if ( i == 5 ) { // x-z
        normal_1 = 1;
        normal_2 = 3;
    } else if ( i == 6 ) { // x-y
        normal_1 = 1;
        normal_2 = 2;
    } else {
        OOFEM_ERROR("Unexpected number for shear stress (must be either 4, 5 or 6).");
        normal_1 = normal_2 = 0;
    }

    if ( ( status->giveTempCrackStatus(normal_1) != pscm_NONE ) && ( status->giveTempCrackStatus(normal_1) != pscm_CLOSED ) ) {
        return false;
    } else if ( ( status->giveTempCrackStatus(normal_2) != pscm_NONE ) && ( status->giveTempCrackStatus(normal_2) != pscm_CLOSED ) ) {
        return false;
    } else {
        return true;
    }
}



double
FCMMaterial :: computeNormalCrackOpening(GaussPoint *gp, int i) {
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    double crackOpening, N;

    if ( i > status->giveNumberOfTempCracks() ) {
        crackOpening = 0;
    } else {
        crackOpening = max(status->giveCharLength(i) * status->giveTempCrackStrain(i), 0.);
        N = this->giveNumberOfCracksInDirection(gp, i);
        crackOpening /= N;
    }

    return crackOpening;
}



double
FCMMaterial :: computeMaxNormalCrackOpening(GaussPoint *gp, int i) {
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    double crackOpening, N;

    if ( i > status->giveNumberOfTempCracks() ) {
        crackOpening = 0;
    } else {
        crackOpening = max(status->giveCharLength(i) * status->giveMaxCrackStrain(i), 0.);
        N = this->giveNumberOfCracksInDirection(gp, i);
        crackOpening /= N;
    }

    return crackOpening;
}



double
FCMMaterial :: computeShearSlipOnCrack(GaussPoint *gp, int icrack) {
    MaterialMode mMode = gp->giveMaterialMode();

    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    // total shear slip
    double slip = 0;

    // number of cracks in the direction of iCrack
    int nCracks = this->giveNumberOfCracksInDirection(gp, icrack);

    // shear directions on icrack plane
    int dir_j, dir_k;

    // cracking shear strains
    double gamma_cr_ij, gamma_cr_ik;

    // factor for redistribution of shear crack strain
    double factor_ij, factor_ik;

    double u_ij, u_ik;

    if ( status->giveTempCrackStatus(icrack) == pscm_NONE ) { // crack not initiated
        return slip;
    }

    if ( icrack > 3 ) {
        OOFEM_ERROR("Unexpected value of index i (4, 5, 6 permitted only)");
    }


    if ( ( mMode == _PlaneStress ) || ( mMode == _PlaneStrain ) ) {
        // get shear strain
        if ( mMode == _PlaneStress ) {
            gamma_cr_ij = status->giveTempCrackStrain(3);
        } else {
            gamma_cr_ij = status->giveTempCrackStrain(6);
        }

        factor_ij = 1.;

        // determine if it necessary to split the cracking strain
        // both cracks exist in 2D

        if ( ( icrack == 2 ) || ( status->giveTempCrackStatus(2) != pscm_NONE ) ) {
            if ( icrack == 1 ) {
                dir_j = 2;
            } else { // icrack == 2
                dir_j = 1;
            }


            // well this is quite unfortunate. The problem is that the shear slip should be redistributed according to the D2 moduli for the individual crack. In reality this is a big problem for the FRC-FCM problem beacuse it would have led to recursive calling (D2 modulus evaluates damage and damage is computed from shear slip etc.

            factor_ij = this->computeShearStiffnessRedistributionFactor(gp, icrack, dir_j);
        }

        slip = factor_ij * fabs(gamma_cr_ij) * status->giveCharLength(icrack) / nCracks;
    } else if ( mMode == _3dMat ) {
        if ( icrack == 1 ) {
            dir_j = 2;
            dir_k = 3;

            gamma_cr_ij = status->giveTempCrackStrain(6);
            gamma_cr_ik = status->giveTempCrackStrain(5);
        } else if ( icrack == 2 ) {
            dir_j = 1;
            dir_k = 3;

            gamma_cr_ij = status->giveTempCrackStrain(6);
            gamma_cr_ik = status->giveTempCrackStrain(4);
        } else { // icrack == 3
            dir_j = 1;
            dir_k = 2;

            gamma_cr_ij = status->giveTempCrackStrain(5);
            gamma_cr_ik = status->giveTempCrackStrain(4);
        }


        factor_ij = 1.;
        factor_ik = 1.;


        if ( ( status->giveTempCrackStatus(dir_j) != pscm_NONE ) || ( status->giveTempCrackStatus(dir_k) != pscm_NONE ) ) {
            if ( status->giveTempCrackStatus(dir_j) != pscm_NONE ) {
                factor_ij = this->computeShearStiffnessRedistributionFactor(gp, icrack, dir_j);
            }

            if ( status->giveTempCrackStatus(dir_k) != pscm_NONE ) {
                factor_ik = this->computeShearStiffnessRedistributionFactor(gp, icrack, dir_k);
            }
        }

        u_ij = factor_ij * fabs(gamma_cr_ij) * status->giveCharLength(icrack) / nCracks;
        u_ik = factor_ik * fabs(gamma_cr_ik) * status->giveCharLength(icrack) / nCracks;

        slip = sqrt( pow(u_ij, 2) + pow(u_ik, 2) );
    } else {
        OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
    }

    return slip;
}


bool
FCMMaterial :: isStrengthExceeded(const FloatMatrix &base, GaussPoint *gp, TimeStep *tStep, int iCrack, double trialStress) {
    if ( trialStress > this->giveTensileStrength(gp) ) {
        return true;
    } else {
        return false;
    }
}




double
FCMMaterial :: computeShearStiffnessRedistributionFactor(GaussPoint *gp, int ithCrackPlane, int jthCrackDirection) {
    double factor_ij;
    double D2_i, D2_j;

    D2_i = this->computeD2ModulusForCrack(gp, ithCrackPlane);
    D2_j = this->computeD2ModulusForCrack(gp, jthCrackDirection);

    factor_ij = D2_j / ( D2_i + D2_j );

    return factor_ij;
}


bool
FCMMaterial :: checkStrengthCriterion(FloatMatrix &newBase, const FloatArray &globalStress, GaussPoint *gp, TimeStep *tStep, int nCrack) {
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    double sigX, sigY, tau, sig2;
    FloatArray trialStress, planeStress, princStress, crackingStrain, shearStrains, newShearStrains;
    FloatMatrix sigmaG2L, princCrackDir, oldBase, princCrackDirExt;

    sigmaG2L = status->giveG2LStressVectorTransformationMtrx();

    trialStress = globalStress;

    // rotate stress to local coordinates
    trialStress.rotatedWith(sigmaG2L, 'n');

    if ( status->giveMaxNumberOfCracks(gp) < 3 ) { // plane stress
        // test material strength but keep crack coordinates
        if ( this->isStrengthExceeded( newBase, gp, tStep, nCrack, trialStress.at(nCrack) ) ) {
            return false;
        } else {
            return true;
        }
    } else if ( nCrack == 3 ) { // 3D but is the last crack
        if ( this->isStrengthExceeded( newBase, gp, tStep, nCrack, trialStress.at(nCrack) ) ) {
            return false;
        } else {
            return true;
        }
    } else { // second crack in 3D case
        // if the stress-strength criterion is violated
        // the crack coordinates have to be changed
        // in-plane shear strain remains unchanged

        // compute second principal stress cheap - if passed, compute eigenvectors and eigenvalues
        sigX = trialStress.at(2); //1st normal stress in crack plane
        sigY = trialStress.at(3); //2nd normal stress in crack plane
        tau = trialStress.at(4); // shear stress in crack plane

        sig2 = ( sigX + sigY ) / 2. + sqrt( ( sigX - sigY ) * ( sigX - sigY ) / 4. + tau * tau );

        if ( this->isStrengthExceeded(newBase, gp, tStep, 2, sig2) ) {
            oldBase = newBase;
            planeStress.resize(3);
            planeStress.zero();

            planeStress.at(1) = trialStress.at(2); //1st normal stress in crack plane
            planeStress.at(2) = trialStress.at(3); //2nd normal stress in crack plane
            planeStress.at(3) = trialStress.at(4); //shear in crack plane


            // compute principal stresses on the already existing crack plane
            this->computePrincipalValDir(princStress, princCrackDir, planeStress, principal_stress);

            // right-hand orientation
            princCrackDir.at(1, 2) = -1. * princCrackDir.at(2, 1);
            princCrackDir.at(2, 2) = princCrackDir.at(1, 1);

            // establish vector of crack directions on the crack plane
            // the first vector is the normal direction, [1, 0, 0]

            princCrackDirExt.resize(3, 3);
            princCrackDirExt.zero();
            princCrackDirExt.at(1, 1) = 1.;

            for ( int i = 1; i <= 2; i++ ) {
                for ( int j = 1; j <= 2; j++ ) {
                    princCrackDirExt.at(i + 1, j + 1) = princCrackDir.at(i, j);
                }
            }

            // new crack base vector in global coordinates
            newBase.beProductOf(oldBase, princCrackDirExt);

            /*
             * // algorithm which gives the same results - rotation of axes in space
             * double theta;
             *
             * if ( princCrackDir.at(1,1) > 0. ) {
             * theta = asin( princCrackDir.at(2,1) );
             * } else {
             * theta = M_PI - asin ( princCrackDir.at(2,1) );
             * }
             *
             * if (theta < 0.) {
             * theta += 2 * M_PI;
             * }
             *
             *
             * FloatArray baseVector2;
             * FloatArray baseVector3;
             *
             * FloatArray newBaseVector2;
             * FloatArray newBaseVector3;
             *
             * baseVector2.resize(3);
             * baseVector3.resize(3);
             *
             * newBaseVector2.resize(3);
             * newBaseVector3.resize(3);
             *
             * for (int i = 1; i <= 3; i++) {
             * baseVector2.at(i) = oldBase.at(i,2);
             * baseVector3.at(i) = oldBase.at(i,3);
             * }
             *
             * FloatMatrix K;
             * FloatMatrix R;
             * FloatMatrix help;
             *
             * K.resize(3,3);
             * K.zero();
             *
             * K.at(1,2) = -1. * oldBase.at(3,1);
             * K.at(1,3) = oldBase.at(2,1);
             * K.at(2,1) = oldBase.at(3,1);
             * K.at(2,3) = -1. * oldBase.at(1,1);
             * K.at(3,1) = -1. * oldBase.at(2,1);
             * K.at(3,2) = oldBase.at(1,1);
             *
             * R.resize(3,3);
             * R.beUnitMatrix();
             *
             * help = K;
             * help.times( sin(theta) );
             *
             * R.add(help);
             *
             * help.beProductOf (K, K);
             * help.times( 1.-cos(theta) );
             *
             * R.add(help);
             *
             * newBaseVector2.beProductOf(R, baseVector2);
             * newBaseVector3.beProductOf(R, baseVector3);
             */


            /*
             * //	the same way - R evaluated symbolically in Matlab
             * //	c = cos, s = sin
             * //	k = axis of rotation
             * R.at(1,1) = 1. + (c-1.) * k2*k2 + (c-1.) * k3*k3;
             * R.at(1,2) = -k3*s - k1*k2 * (c-1.);
             * R.at(1,3) =  k2*s - k1*k3 * (c-1.);
             *
             * R.at(2,1) = k3*s - k1*k2 * (c-1.);
             * R.at(2,2) = 1. + (c-1.) * k1*k1 + (c-1.) * k3*k3;
             * R.at(2,3) = -k1*s - k2*k3 * (c-1.);
             *
             * R.at(3,1) = -k2*s - k1*k3 * (c-1.);
             * R.at(3,2) = k1*s - k2*k3 * (c-1.);
             * R.at(3,3) = 1. + (c-1.) * k1*k1 + (c-1.) * k2*k2;
             */


            // modify shear strains

            crackingStrain = status->giveCrackStrainVector();

            shearStrains.resize(2);
            shearStrains.at(1) = crackingStrain.at(6); //aka x-y
            shearStrains.at(2) = crackingStrain.at(5); //aka x-z

            // gamma_xy_l   [ c   s ]    gamma_xy_g
            //            = [       ] *
            // gamma_xz_l   [ -s  c ]    gamma_xz_g

            newShearStrains.beTProductOf(princCrackDir, shearStrains);

            crackingStrain.at(5) =  newShearStrains.at(2); // aka x-z
            crackingStrain.at(6) =  newShearStrains.at(1); // aka x-y

            status->setCrackStrainVector(crackingStrain);
            status->setTempCrackStrainVector(crackingStrain);

            // modify maximum shear strains
            // get max shear strains in the original configuration
            shearStrains.at(1) = status->giveMaxCrackStrain(6); //aka x-y
            shearStrains.at(2) = status->giveMaxCrackStrain(5); //aka x-z

            // gamma_xy_l   [ c   s ]    gamma_xy_g
            //            = [       ] *
            // gamma_xz_l   [ -s  c ]    gamma_xz_g

            newShearStrains.beTProductOf(princCrackDir, shearStrains);

            status->setMaxCrackStrain( 5, max( newShearStrains.at(2), crackingStrain.at(5) ) ); // aka x-z
            status->setMaxCrackStrain( 6, max( newShearStrains.at(1), crackingStrain.at(6) ) ); // aka x-y

            status->setTempMaxCrackStrain( 5, max( newShearStrains.at(2), crackingStrain.at(5) ) ); // aka x-z
            status->setTempMaxCrackStrain( 6, max( newShearStrains.at(1), crackingStrain.at(6) ) ); // aka x-y

            return false;
        } else { // strength not reached
            return true;
        }
    } // number of crack in given stress-state
}


double
FCMMaterial :: giveCharacteristicElementLength(GaussPoint *gp, const FloatArray &crackPlaneNormal)
{
    return gp->giveElement()->giveCharacteristicLength(crackPlaneNormal);
}


void
FCMMaterial :: updateCrackStatus(GaussPoint *gp)
{
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    int nMaxCracks = status->giveMaxNumberOfCracks(gp);

    FloatArray crackStrain, maxCrackStrain;
    crackStrain = status->giveTempCrackStrainVector();
    maxCrackStrain = status->giveMaxCrackStrainVector();


    // loop over NORMAL components of the crack strain vector
    // statuses are only for normal components
    for ( int i = 1; i <= nMaxCracks; i++ ) {
        // changes in status only in previously cracked (or at least initiated in this step) material
        if ( status->giveTempCrackStatus(i) != pscm_NONE ) {
            // temp crack bigger or equal than max strain so far
            if ( crackStrain.at(i) >= maxCrackStrain.at(i) ) {
                status->setTempMaxCrackStrain( i, crackStrain.at(i) );

                if ( status->giveCrackStatus(i) != pscm_NONE ) { //already existing crack in prev step
                    status->setTempCrackStatus(i, pscm_SOFTENING);
                } else if ( ( fabs( crackStrain.at(i) ) <= 1.e-18 ) && ( status->giveCrackStatus(i) == pscm_NONE ) ) { // just initiated crack with zero strain
                    if ( i == nMaxCracks ) {
                        status->setTempCrackStatus(i, pscm_NONE);
                    } else if ( nMaxCracks > i ) { // curent crack is not the last possible one
                        if ( status->giveTempCrackStatus(i + 1) == pscm_NONE ) { // the next crack is not initiated
                            status->setTempCrackStatus(i, pscm_NONE); // reset the current one to NONE
                        } else {
                            status->setTempCrackStatus(i, pscm_JUST_INIT); // else keep as just initiated
                        }
                    } else {
                        status->setTempCrackStatus(i, pscm_JUST_INIT);
                    }
                } else {
                    status->setTempCrackStatus(i, pscm_JUST_INIT);
                }

                // closed previously
            } else if ( crackStrain.at(i) <= 0. ) {
                status->setTempCrackStatus(i, pscm_CLOSED);

                if ( status->giveMaxCrackStrain(i) == 0. ) { //not existing crack in previous step that wants to close?
                    status->setTempMaxCrackStrain(i, 0.);
                    //	  status->setTempCrackStatus(i, pscm_NONE);
                    // CLOSED status should have similar behavior and it is safer. The status can be changed within one iteration loop on the local scale
                    // the crack can be changed from CLOSED to NONE during overall updating
                    //	  status->setTempCrackStatus(i, pscm_CLOSED);
                }

                // unloading--reloading
            } else if ( crackStrain.at(i) < maxCrackStrain.at(i) ) {
                status->setTempCrackStatus(i, pscm_UNLO_RELO);
            } else {
                OOFEM_ERROR("Unexpected value of cracking strain");
            }
        }
    }

    // SHEAR components
    for ( int i = nMaxCracks + 1; i <=  crackStrain.giveSize(); i++ ) {
        if ( fabs( crackStrain.at(i) ) > maxCrackStrain.at(i) ) {
            status->setTempMaxCrackStrain( i, fabs( crackStrain.at(i) ) );
        }
    }
}





void
FCMMaterial :: giveLocalCrackedStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode, GaussPoint *gp,
                                               TimeStep *tStep)
{
    int dim, j;
    FloatMatrix Dcr;
    IntArray indx;
    StructuralMaterial :: giveVoigtSymVectorMask( indx, gp->giveMaterialMode() );

    dim = indx.giveSize();

    Dcr.resize(dim, dim);

    for ( int i = 1; i <= dim; i++ ) {
        j = indx.at(i);

        if ( j <= 3 ) {
            Dcr.at(i, i) = this->giveCrackingModulus(rMode, gp, j);
        } else {
            Dcr.at(i, i) = this->computeTotalD2Modulus(gp, j);
        }
    }

    answer = Dcr;
}





void
FCMMaterial :: giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                           MatResponseMode rMode,
                                           GaussPoint *gp,
                                           TimeStep *tStep)
//
// returns effective material stiffness matrix in full form
// for gp stress strain mode
//
{
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );
    StructuralMaterial *lMat = static_cast< StructuralMaterial * >( this->giveLinearElasticMaterial() );

    MaterialMode mMode = gp->giveMaterialMode();

    int numberOfActiveCracks, nMaxCracks;
    double overallElasticStiffness;

    FloatMatrix D, De, DeHelp, Dcr, DcrHelp, DcrHelp2, inv, stiffnessL2G;

    numberOfActiveCracks = status->giveNumberOfTempCracks();

    // ELASTIC MATRIX
    if ( ( rMode == ElasticStiffness ) || ( numberOfActiveCracks == 0 ) ) {
        lMat->giveStiffnessMatrix(D, rMode, gp, tStep);

        overallElasticStiffness = this->computeOverallElasticStiffness();
        if ( overallElasticStiffness != ( this->give('E', gp) ) ) {
            D.times( overallElasticStiffness / ( this->give('E', gp) ) );
        }

        answer = D;
        return;
    }

    // SECANT OR TANGENT MATRIX - first in local direction
    nMaxCracks = status->giveMaxNumberOfCracks(gp);
    lMat->giveStiffnessMatrix(De, rMode, gp, tStep);

    overallElasticStiffness = this->computeOverallElasticStiffness();
    if ( overallElasticStiffness != ( this->give('E', gp) ) ) {
        De.times( overallElasticStiffness / ( this->give('E', gp) ) );
    }

    // extract 1x1 / 2x2 / 3x3 submatrix
    De.resizeWithData(nMaxCracks, nMaxCracks);

    // shrink to square matrix number of cracks x number of cracks
    DeHelp = De;
    DeHelp.resizeWithData(numberOfActiveCracks, numberOfActiveCracks);

    Dcr.resize(numberOfActiveCracks, numberOfActiveCracks);
    Dcr.zero();

    for ( int i = 1; i <= numberOfActiveCracks; i++ ) {
        Dcr.at(i, i) = this->giveCrackingModulus(rMode, gp, i);
    }

    Dcr.add(DeHelp);
    inv.beInverseOf(Dcr);
    inv.resizeWithData(nMaxCracks, nMaxCracks);


    DcrHelp.beProductOf(De, inv); // De (De + Dcr)^-1
    DcrHelp2.beProductOf(DcrHelp, De); // De (De + Dcr)^-1 De

    D.zero();
    D.resize(nMaxCracks, nMaxCracks);
    D.add(De);
    D.subtract(DcrHelp2); // De - De (De + Dcr)^-1 De


    // resize add shear moduli on diagonal
    switch ( mMode ) {
    case _PlaneStress:
        D.resizeWithData(3, 3);
        D.at(3, 3) = this->computeEffectiveShearModulus(gp, 6);

        break;

    case _3dMat:
        D.resizeWithData(6, 6);
        for ( int i = 4; i <= 6; i++ ) {
            D.at(i, i) = this->computeEffectiveShearModulus(gp, i);
        }
        break;

    case _PlaneStrain:
        D.resizeWithData(6, 6);
        D.at(6, 6) = this->computeEffectiveShearModulus(gp, 6);
        break;

    default:
        OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
    }

    // transform stiffnes to global c.s
    // to transform stiffness from LOCAL to GLOBAL coordinate system
    // the transformation matrix is the same as for strain transformation
    // from GLOBAL to LOCAL coordinate system
    stiffnessL2G = status->giveG2LStrainVectorTransformationMtrx();
    D.rotatedWith(stiffnessL2G, 'n');


    switch ( mMode ) {
    case _PlaneStress:
        answer = D;
        return;

        break;
    case _3dMat:
        answer = D;
        return;

        break;
    case _PlaneStrain:
        StructuralMaterial :: giveReducedSymMatrixForm( answer, D, gp->giveMaterialMode() );
        return;

        break;
    default:
        OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
        return;
    }
}



double
FCMMaterial :: computeTotalD2Modulus(GaussPoint *gp, int shearDirection)
{
    double D2 = 0.;
    int crackA, crackB;
    double D2_1, D2_2;

    if ( this->isIntactForShear(gp, shearDirection) ) {
        D2 = this->computeOverallElasticStiffness() * fcm_BIGNUMBER;
    } else {
        if ( shearDirection == 4 ) {
            crackA = 2;
            crackB = 3;
        } else if ( shearDirection == 5 ) {
            crackA = 1;
            crackB = 3;
        }  else if ( shearDirection == 6 ) {
            crackA = 1;
            crackB = 2;
        } else {
            OOFEM_ERROR("Unexpected value of index i (4, 5, 6 permitted only)");
            crackA = crackB = 0;
        }

        if ( ( this->isIntact(gp, crackA) ) || ( this->isIntact(gp, crackB) ) ) {
            if ( this->isIntact(gp, crackA) ) {
                D2 = this->computeD2ModulusForCrack(gp, crackB);
            } else {
                D2 = this->computeD2ModulusForCrack(gp, crackA);
            }
        } else {
            D2_1 = this->computeD2ModulusForCrack(gp, crackA);
            D2_2 = this->computeD2ModulusForCrack(gp, crackB);

            if ( multipleCrackShear ) {
                // serial coupling of stiffnesses 1/D = 1/D1 + 1/D2
                D2 = D2_1 * D2_2 / ( D2_1 + D2_2 );
            } else {
                D2 = min(D2_1, D2_2);
            }
        }
    }

    return D2;
}

IRResultType
FCMMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    this->nAllowedCracks = 3;
    IR_GIVE_OPTIONAL_FIELD(ir, nAllowedCracks, _IFT_FCM_nAllowedCracks);


    this->crackSpacing = -1.;
    if  ( ir->hasField(_IFT_FCM_crackSpacing) ) {
        IR_GIVE_FIELD(ir, crackSpacing, _IFT_FCM_crackSpacing);
    }

    this->multipleCrackShear = false;
    if  ( ir->hasField(_IFT_FCM_multipleCrackShear) ) {
        this->multipleCrackShear = true;
    }

    int ecsm = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, ecsm, _IFT_FCM_ecsm);
    switch ( ecsm ) {
    case 1: ecsMethod = ECSM_SquareRootOfArea;
        break;
    case 2: ecsMethod = ECSM_ProjectionCentered;
        break;
    case 3: ecsMethod = ECSM_Oliver1;
        break;
    case 4: ecsMethod = ECSM_Oliver1modified;
        break;
    default: ecsMethod = ECSM_Projection;
    }


    return IRRT_OK;
}


double
FCMMaterial :: give(int aProperty, GaussPoint *gp)
{
    return linearElasticMaterial->give(aProperty, gp);
}



double
FCMMaterial :: giveCrackSpacing(void)
{
    return crackSpacing;
}


double
FCMMaterial :: giveNumberOfCracksInDirection(GaussPoint *gp, int iCrack)
{
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );
    double spacing, L, N;

    L = status->giveCharLength(iCrack);
    spacing = this->giveCrackSpacing();

    if ( ( spacing > L ) || ( spacing < 0. ) ) {
        N = 1;
    } else {
        N = ( L / spacing );
    }

    return N;
}


double
FCMMaterial :: giveNumberOfCracksForShearDirection(GaussPoint *gp, int i)
{
    double N;
    int dir_1, dir_2;

    if ( i == 4 ) {
        dir_1 = 2;
        dir_2 = 3;
    } else if ( i == 5 ) {
        dir_1 = 1;
        dir_2 = 3;
    }  else if ( i == 6 ) {
        dir_1 = 1;
        dir_2 = 2;
    } else {
        OOFEM_ERROR("Unexpected value of index i (4, 5, 6 permitted only)");
        dir_1 = dir_2 = 0;
    }

    N = max( this->giveNumberOfCracksInDirection(gp, dir_1), this->giveNumberOfCracksInDirection(gp, dir_2) );

    return N;
}


int
FCMMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );
    double width;
    int index;


    // first/dominant crack vector
    if ( type == IST_CrackVector ) {
        answer.resize(3);
        answer.zero();

        // MAX CRACK
        width = 0.;
        index = 1;
        for ( int i = 1; i <= status->giveNumberOfCracks(); i++ ) {
            if ( status->giveCharLength(i) * status->giveCrackStrain(i) / this->giveNumberOfCracksInDirection(gp, i)  > width ) {
                width = status->giveCharLength(i) * status->giveCrackStrain(i) / this->giveNumberOfCracksInDirection(gp, i);
                index = i;
            }
        }

        for ( int i = 1; i <= status->giveCrackDirs().giveNumberOfRows(); i++ ) {
            answer.at(i) = status->giveCrackDirs().at(i, index);
        }

        return 1;
    } else if ( type == IST_2ndCrackVector ) {
        answer.resize(3);
        answer.zero();
        index = 2;

        if ( status->giveNumberOfCracks() <= 1 ) {
            return 1;
        } else if ( status->giveNumberOfCracks() == 2 ) {
            if ( status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) {
                index = 1;
            }
        } else { // 3 cracks
            if ( status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) { // #2 is biggest
                if ( status->giveCharLength(3) * status->giveCrackStrain(3) / this->giveNumberOfCracksInDirection(gp, 3) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) {
                    index = 3;
                } else {
                    index = 1;
                }
            } else { // #1 is biggest
                if ( status->giveCharLength(3) * status->giveCrackStrain(3) / this->giveNumberOfCracksInDirection(gp, 3) > status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) ) {
                    index = 3;
                } else {
                    index = 2;
                }
            }
        }

        for ( int i = 1; i <= status->giveCrackDirs().giveNumberOfRows(); i++ ) {
            answer.at(i) = status->giveCrackDirs().at(i, index);
        }

        return 1;
    } else if ( type == IST_3rdCrackVector ) {
        answer.resize(3);
        answer.zero();
        index = 3;

        if ( status->giveNumberOfCracks() <= 2 ) {
            return 1;
        } else { // 3 cracks
            if ( status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) { // #2 is biggest
                if ( status->giveCharLength(3) * status->giveCrackStrain(3) / this->giveNumberOfCracksInDirection(gp, 3) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) {
                    index = 1;
                } else {
                    index = 3;
                }
            } else { // #1 is biggest
                if ( status->giveCharLength(3) * status->giveCrackStrain(3) / this->giveNumberOfCracksInDirection(gp, 3) > status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) ) {
                    index = 2;
                } else {
                    index = 3;
                }
            }
        }

        for ( int i = 1; i <= status->giveCrackDirs().giveNumberOfRows(); i++ ) {
            answer.at(i) = status->giveCrackDirs().at(i, index);
        }

        return 1;


        // width of a first/dominant crack
    } else if ( type == IST_CrackWidth ) {
        answer.resize(1);
        answer.zero();

        // MAX WIDTH
        width = 0.;
        for ( int i = 1; i <= status->giveNumberOfCracks(); i++ ) {
            width = max( width, status->giveCharLength(i) * status->giveCrackStrain(i) / this->giveNumberOfCracksInDirection(gp, i) );
        }
        answer.at(1) = width;
        return 1;
    } else if ( type == IST_2ndCrackWidth ) {
        answer.resize(1);
        answer.zero();
        index = 2;
        if ( status->giveNumberOfCracks() <= 1 ) {
            return 1;
        } else if ( status->giveNumberOfCracks() == 2 ) {
            if ( status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2)  > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) {
                index = 1;
            }
        } else { // 3 cracks
            if ( status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) { // #2 is biggest
                if ( status->giveCharLength(3) * status->giveCrackStrain(3) / this->giveNumberOfCracksInDirection(gp, 3) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) {
                    index = 3;
                } else {
                    index = 1;
                }
            } else { // #1 is biggest
                if ( status->giveCharLength(3) * status->giveCrackStrain(3) / this->giveNumberOfCracksInDirection(gp, 3) > status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) ) {
                    index = 3;
                } else {
                    index = 2;
                }
            }
        }

        width = status->giveCharLength(index) * status->giveCrackStrain(index) / this->giveNumberOfCracksInDirection(gp, index);

        answer.at(1) = width;
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
    } else if ( type == IST_CrackStrainTensor ) {
        FloatArray crackStrain = status->giveCrackStrainVector();
        FloatMatrix epsL2G = status->giveL2GStrainVectorTransformationMtrx();
        // from local to global
        crackStrain.rotatedWith(epsL2G, 'n');

        StructuralMaterial :: giveFullSymVectorForm( answer, crackStrain, gp->giveMaterialMode() );

        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

void
FCMMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                             MatResponseMode mode,
                                             GaussPoint *gp,
                                             TimeStep *tStep)
{
    this->giveMaterialStiffnessMatrix(answer, mode, gp, tStep);
}


void
FCMMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode mode,
                                        GaussPoint *gp,
                                        TimeStep *tStep)

{
    this->giveMaterialStiffnessMatrix(answer, mode, gp, tStep);
}


void
FCMMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode mode,
                                        GaussPoint *gp,
                                        TimeStep *tStep)

{
    this->giveMaterialStiffnessMatrix(answer, mode, gp, tStep);
}


FCMMaterialStatus :: FCMMaterialStatus(int n, Domain *d, GaussPoint *gp) :
    StructuralMaterialStatus(n, d, gp),
    crackStatuses(), tempCrackStatuses(),
    maxCrackStrains(), tempMaxCrackStrains(),
    crackStrainVector(), tempCrackStrainVector(),
    crackDirs(),
    charLengths(),
    transMatrix_G2Lstress(), transMatrix_G2Lstrain(),
    transMatrix_L2Gstress(), transMatrix_L2Gstrain()
{
    // resize in constructor according to stress-state
    this->nMaxCracks = 0;
    this->nMaxCracks = giveMaxNumberOfCracks(gp);

    crackStatuses.resize(this->nMaxCracks);
    crackStatuses.zero();
    tempCrackStatuses = crackStatuses;

    charLengths.resize(this->nMaxCracks);
    charLengths.zero();

    crackDirs.resize(this->nMaxCracks, this->nMaxCracks);
    crackDirs.zero();
    for ( int i = 1; i <= this->nMaxCracks; i++ ) {
        crackDirs.at(i, i) = 1.0;
    }


    if ( this->nMaxCracks == 2 ) { //plane stress
        maxCrackStrains.resize(3);
        maxCrackStrains.zero();
        tempMaxCrackStrains = maxCrackStrains;

        crackStrainVector = maxCrackStrains;
        tempCrackStrainVector = maxCrackStrains;


        transMatrix_G2Lstress.resize(3, 3);
        transMatrix_G2Lstress.zero();
        transMatrix_L2Gstrain = transMatrix_L2Gstress = transMatrix_G2Lstrain = transMatrix_G2Lstress;
    } else {
        maxCrackStrains.resize(6);
        maxCrackStrains.zero();
        tempMaxCrackStrains = maxCrackStrains;

        crackStrainVector = maxCrackStrains;
        tempCrackStrainVector = maxCrackStrains;


        transMatrix_G2Lstress.resize(6, 6);
        transMatrix_G2Lstress.zero();
        transMatrix_L2Gstrain = transMatrix_L2Gstress = transMatrix_G2Lstrain = transMatrix_G2Lstress;
    }
}


FCMMaterialStatus :: ~FCMMaterialStatus()
{ }



void
FCMMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    int i;
    char s [ 11 ];

    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( this->giveNumberOfCracks() > 0 ) {
        for ( i = 1; i <= crackDirs.giveNumberOfColumns(); i++ ) {
            switch ( crackStatuses.at(i) ) {
            case pscm_NONE:
                strcpy(s, "NONE");
                break;
            case pscm_JUST_INIT:
                strcpy(s, "JUST_INIT");
                break;
            case pscm_SOFTENING:
                strcpy(s, "SOFTENING");
                break;
            case pscm_UNLO_RELO:
                strcpy(s, "UNLO_RELO");
                break;
            case pscm_CLOSED:
                strcpy(s, "CLOSED");
                break;
            default:
                strcpy(s, "UNKNOWN");
                break;
            }

            fprintf(file, "crack %d {status %s, crackplane is normal to { ", i, s);

            for ( int j = 1; j <= crackDirs.giveNumberOfRows(); j++ ) {
                fprintf( file, "%f ", crackDirs.at(j, i) );
            }

            fprintf(file, "}}");
            ;
        }
    }

    fprintf(file, " }\n");
}

int
FCMMaterialStatus :: giveNumberOfCracks() const
//
// return number of existing cracks
//
{
    int answer = 0;

    for ( int i = 1; i <= crackStatuses.giveSize(); i++ ) {
        if ( crackStatuses.at(i) != pscm_NONE  ) {
            answer++;
        }
    }

    return answer;
}


int
FCMMaterialStatus :: giveNumberOfTempCracks() const
//
// return number of existing temp cracks
//
{
    int answer = 0;

    for ( int i = 1; i <= tempCrackStatuses.giveSize(); i++ ) {
        if ( tempCrackStatuses.at(i) != pscm_NONE  ) {
            answer++;
        }
    }

    return answer;
}

int
FCMMaterialStatus :: giveMaxNumberOfCracks(GaussPoint *gp)
//
// return number of maximum allowable cracks
//
{
    int nCr = 0;
    IntArray indx;

    if ( this->nMaxCracks == 0 ) {
        StructuralMaterial :: giveVoigtSymVectorMask( indx, gp->giveMaterialMode() );

        for ( int i = 1; i <= 3; i++ ) {
            if ( indx.contains(i) ) {
                nCr++;
            }
        }

        this->nMaxCracks = nCr;
    }

    return this->nMaxCracks;
}




void
FCMMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
//
{
    StructuralMaterialStatus :: initTempStatus();

    tempCrackStatuses = crackStatuses;
    tempMaxCrackStrains = maxCrackStrains;
    tempCrackStrainVector = crackStrainVector;
}



void
FCMMaterialStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables correspond to newly reched equilibrium.
//
{
    StructuralMaterialStatus :: updateYourself(tStep);


    maxCrackStrains = tempMaxCrackStrains;
    crackStrainVector = tempCrackStrainVector;

    for ( int i = 1; i <= crackStrainVector.giveSize(); i++ ) {
        if ( crackStrainVector.at(i) < 0. ) {
            crackStrainVector.at(i) = 0.;
        }
    }

    // updating of statuses has to be done more carefully.
    // consider a crack which does not exist in the previous step and ends as closed in the end of this step
    // this crack is naturally treated as "NONE" in the following steps

    for ( int i = 1; i <= crackStatuses.giveSize(); i++ ) {
        if ( ( tempCrackStatuses.at(i) == pscm_CLOSED ) && ( crackStatuses.at(i) == pscm_NONE ) ) {
            // no other crack so this one can be set as non-existing
            if ( i + 1 > nMaxCracks ) {
                crackStatuses.at(i) = pscm_NONE;


                // be sure that in the second and third crack does not exist, if it does we have to copy CLOSED status
            } else if ( tempCrackStatuses.at(i + 1) != pscm_NONE ) {
                crackStatuses.at(i) = tempCrackStatuses.at(i);
            } else {
                crackStatuses.at(i) = pscm_NONE;
            }
        } else {
            crackStatuses.at(i) = tempCrackStatuses.at(i);
        }
    }
}



contextIOResultType
FCMMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
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

    if ( ( iores = charLengths.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackStrainVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackStrainVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_G2Lstrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_G2Lstress.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_L2Gstrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_L2Gstress.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

contextIOResultType
FCMMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
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

    if ( ( iores = charLengths.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackStrainVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_G2Lstrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_G2Lstress.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_L2Gstrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_L2Gstress.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK; // return succes
}
} // end namespace oofem
