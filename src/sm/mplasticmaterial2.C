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

#include "mplasticmaterial2.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"

#include "structuralcrosssection.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
#define YIELD_TOL 1.e-6
#define YIELD_TOL_SECONDARY 1.e-4
#define RES_TOL   1.e-8
#define PLASTIC_MATERIAL_MAX_ITERATIONS 120

MPlasticMaterial2 :: MPlasticMaterial2(int n, Domain *d)  : StructuralMaterial(n, d)
    //
    // constructor
    //
{
    linearElasticMaterial = NULL;
    rmType = mpm_ClosestPoint;
    plType = associatedPT;
    iterativeUpdateOfActiveConds = false; // safe
}



MPlasticMaterial2 :: ~MPlasticMaterial2()
//
// destructor
//
{
    delete linearElasticMaterial;
}


int
MPlasticMaterial2 :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    if ( ( mode == _3dMat ) ||
        ( mode == _1dMat ) ||
        //<RESTRICTED_SECTION>
        ( mode == _PlaneStress )  ||
        //</RESTRICTED_SECTION>
        ( mode == _PlaneStrain )  ||
        ( mode == _2dPlateLayer ) ||
        ( mode == _2dBeamLayer )  ||
        ( mode == _3dShellLayer ) ||
        ( mode == _1dFiber ) ) {
        return 1;
    }

    return 0;
}

MaterialStatus *
MPlasticMaterial2 :: CreateStatus(GaussPoint *gp) const
/*
 * creates new  material status  corresponding to this class
 */
{
    MPlasticMaterial2Status *status;

    status = new MPlasticMaterial2Status(1, this->giveDomain(), gp);
    return status;
}

void
MPlasticMaterial2 :: giveRealStressVector(FloatArray &answer,
                                          MatResponseForm form,
                                          GaussPoint *gp,
                                          const FloatArray &totalStrain,
                                          TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
// completely formulated in Reduced stress-strain space
{
    FloatArray strainSpaceHardeningVariables;
    FloatArray fullStressVector;
    FloatArray strainIncrement, elasticStrainVectorR;
    FloatArray strainVectorR, plasticStrainVectorR;
    FloatArray helpVec, helpVec2;
    int i;
    FloatMatrix elasticModuli, hardeningModuli, consistentModuli;
    FloatMatrix elasticModuliInverse, hardeningModuliInverse;
    FloatMatrix helpMtrx, helpMtrx2;
    IntArray activeConditionMap(this->nsurf);
    FloatArray gamma;

    MPlasticMaterial2Status *status = ( MPlasticMaterial2Status * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(strainVectorR, gp, totalStrain,
                                                atTime, VM_Total);

    /*
     * stress return algorithm
     */

    //if (0) this->cuttingPlaneReturn (fullStressVector, activeConditionMap, gamma, form, gp, totalStrain, atTime);
    //else this->closestPointReturn (fullStressVector, activeConditionMap, gamma, form, gp, totalStrain, atTime);
    if ( rmType == mpm_ClosestPoint ) {
        this->closestPointReturn(fullStressVector, activeConditionMap, gamma, form, gp, strainVectorR,
                                 plasticStrainVectorR, strainSpaceHardeningVariables, atTime);
    } else {
        this->cuttingPlaneReturn(fullStressVector, activeConditionMap, gamma, form, gp, strainVectorR,
                                 plasticStrainVectorR, strainSpaceHardeningVariables, atTime);
    }


    status->letTempStrainVectorBe(totalStrain);
    crossSection->giveReducedCharacteristicVector(helpVec, gp, fullStressVector);
    status->letTempStressVectorBe(helpVec);

    status->letTempPlasticStrainVectorBe(plasticStrainVectorR);
    status->letTempStrainSpaceHardeningVarsVectorBe(strainSpaceHardeningVariables);

    status->setTempGamma(gamma);
    status->setTempActiveConditionMap(activeConditionMap);


    // update state flag
    int newState, state = status->giveStateFlag();
    bool yieldFlag = false;
    for ( i = 1; i <= nsurf; i++ ) {
        if ( gamma.at(i) > 0. ) {
            yieldFlag = true;
        }
    }

    if ( yieldFlag ) {
        newState = MPlasticMaterial2Status :: PM_Yielding;         // test if plastic loading occur
    }
    // no plastic loading - check for unloading
    else if ( ( state == MPlasticMaterial2Status :: PM_Yielding ) || ( state == MPlasticMaterial2Status :: PM_Unloading ) ) {
        newState = MPlasticMaterial2Status :: PM_Unloading;
    } else {
        newState = MPlasticMaterial2Status :: PM_Elastic;
    }

    status->letTempStateFlagBe(newState);

    if ( form == FullForm ) {
        answer = fullStressVector;
        return;
    } else {
        crossSection->giveReducedCharacteristicVector(answer, gp, fullStressVector);
        return;
    }
}


void
MPlasticMaterial2 :: closestPointReturn(FloatArray &answer,
                                        IntArray &activeConditionMap,
                                        FloatArray &gamma,
                                        MatResponseForm form,
                                        GaussPoint *gp,
                                        const FloatArray &totalStrain,
                                        FloatArray &plasticStrainVectorR,
                                        FloatArray &strainSpaceHardeningVariables,
                                        TimeStep *atTime)
{
    FloatArray fullStressVector;
    FloatArray strainIncrement, elasticStrainVectorR;
    FloatArray fullStressSpaceHardeningVars, residualVectorR, gradientVectorR;
    ;
    FloatArray helpVector, helpVector2;
    FloatArray dgamma, dgammared, tempGamma, dkappa;
    FloatMatrix elasticModuli, hardeningModuli, consistentModuli;
    FloatMatrix elasticModuliInverse, hardeningModuliInverse;
    FloatMatrix helpMtrx, helpMtrx2;
    FloatMatrix ks, kl, lmat, rmat, gradientMatrix;
    FloatMatrix gmat;
    IntArray initialConditionMap;
    std :: vector< FloatArray >yieldGradSigVec(this->nsurf), loadGradSigVec(this->nsurf), * yieldGradSigVecPtr, * loadGradSigVecPtr;
    std :: vector< FloatArray >yieldGradKVec(this->nsurf), loadGradKVec(this->nsurf), * yieldGradKVecPtr, * loadGradKVecPtr;
    FloatArray rhs;
    double yieldValue;
    int nIterations = 0;
    int i, j, strSize, totSize;
    int elastic, restart, actSurf, indx;
    bool yieldConsistency, init = true;
    bool hasHardening = this->hasHardening();
#ifdef DEBUG
    bool debug = false;
#endif

    this->clearPopulationSet();

    if ( this->plType == associatedPT ) {
        yieldGradSigVecPtr = loadGradSigVecPtr = & yieldGradSigVec;
        yieldGradKVecPtr = loadGradKVecPtr = & yieldGradKVec;
    } else {
        yieldGradSigVecPtr = & yieldGradSigVec;
        loadGradSigVecPtr  = & loadGradSigVec;
        yieldGradKVecPtr = & yieldGradKVec;
        loadGradKVecPtr  = & loadGradKVec;
    }

    MPlasticMaterial2Status *status = ( MPlasticMaterial2Status * ) this->giveStatus(gp);

huhu:

    status->givePlasticStrainVector(plasticStrainVectorR);
    status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);

    dgamma.resize(nsurf);
    dgamma.zero();
    gamma.resize(nsurf);
    gamma.zero();

    strSize = totalStrain.giveSize(); // size of reducedStrain Vector
    totSize = strSize + strainSpaceHardeningVariables.giveSize();

    // compute elastic moduli and its inverse
    this->computeReducedElasticModuli(elasticModuli, gp, atTime);
    elasticModuliInverse.beInverseOf(elasticModuli);
    // delete elasticModuli;

    //
    // Elastic Predictor
    //

    elasticStrainVectorR = totalStrain;
    elasticStrainVectorR.subtract(plasticStrainVectorR);
    // stress vector in full form due to computational convinience
    //if (fullStressVector) delete fullStressVector;
    this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
    //
    // check for plastic process
    //
    elastic = 1;
    actSurf = 0;
    activeConditionMap.zero();
    if ( init ) {
        initialConditionMap.resize(nsurf);
        initialConditionMap.zero();
    }

    for ( i = 1; i <= this->nsurf; i++ ) {
        yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, strainSpaceHardeningVariables);
        if ( yieldValue > YIELD_TOL ) {
            elastic = 0;
            if ( init ) {
                initialConditionMap.at(i) = 1;
            }

            if ( actSurf < this->giveMaxNumberOfActiveYieldConds(gp) ) {
                actSurf += 1;
                activeConditionMap.at(i) = actSurf;
            }
        }
    }

    init = false;

    if ( !elastic ) {
        /*
         * plastic multisurface closest-point corrector
         */
        this->addNewPopulation(activeConditionMap);

#ifdef DEBUG
        if ( debug ) {
            printf("New activeConditionMap:");
            activeConditionMap.printYourself();
            printf("\n");
        }

#endif

        do {
            do { // restart loop
                 // compute gradients
                for ( i = 1; i <= nsurf; i++ ) {
                    this->computeReducedStressGradientVector(yieldGradSigVec [ i - 1 ], yieldFunction, i, gp, fullStressVector,
                                                             strainSpaceHardeningVariables);
                    if ( hasHardening ) {
                        this->computeKGradientVector(yieldGradKVec [ i - 1 ], yieldFunction, i, gp, fullStressVector,
                                                     strainSpaceHardeningVariables);
                    }
                }

                if ( this->plType == nonassociatedPT ) {
                    for ( i = 1; i <= nsurf; i++ ) {
                        this->computeReducedStressGradientVector(loadGradSigVec [ i - 1 ], loadFunction, i, gp, fullStressVector,
                                                                 strainSpaceHardeningVariables);
                        if ( hasHardening ) {
                            this->computeKGradientVector(loadGradKVec [ i - 1 ], loadFunction, i, gp, fullStressVector,
                                                         strainSpaceHardeningVariables);
                        }
                    }
                }


                //this->computeGradientVector (gradientVectorR, gp, &fullStressVector, fullStressSpaceHardeningVars);
                this->computeResidualVector(residualVectorR, gp, gamma, activeConditionMap, plasticStrainVectorR,
                                            strainSpaceHardeningVariables, * loadGradSigVecPtr);
#ifdef DEBUG
                if ( debug ) {
                    printf("Convergence[actSet: ");
                    for ( i = 1; i <= nsurf; i++ ) {
                        printf( "%d ", activeConditionMap.at(i) );
                    }

                    printf("] yield_val (gamma) ");
                }

#endif

                // check convergence
                yieldConsistency = true;
                for ( i = 1; i <= nsurf; i++ ) {
                    if ( activeConditionMap.at(i) ) {
                        yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, strainSpaceHardeningVariables);
                        if ( fabs(yieldValue) >= YIELD_TOL ) {
                            yieldConsistency = false;
                        }

#ifdef DEBUG
                        if ( debug ) {
                            printf( "%e(%e) ", yieldValue, gamma.at(i) );
                        }

#endif
                    }
                }

#ifdef DEBUG
                if ( debug ) {
                    printf("\n");
                }

#endif

                //if (yieldConsistency && (sqrt(dotProduct(residualVectorR, residualVectorR, residualVectorR.giveSize())) < RES_TOL)) {
                if ( yieldConsistency ) {
                    answer = fullStressVector;
                    //printf (" (%d iterations)", nIterations);


                    restart = 0;
                    // check for active conditions
#ifdef DEBUG
                    if ( debug ) {
                        printf("Consistency: ");
                    }

#endif
                    for ( i = 1; i <= this->nsurf; i++ ) {
                        yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, strainSpaceHardeningVariables);

#ifdef DEBUG
                        if ( debug ) {
                            printf("%e ", yieldValue);
                        }

#endif

                        if ( ( gamma.at(i) < 0.0 ) || ( ( activeConditionMap.at(i) == 0 ) && ( yieldValue > YIELD_TOL_SECONDARY ) ) ) {
                            restart = 1;
                            //printf ("*");
                        }
                    }

#ifdef DEBUG
                    if ( debug ) {
                        printf("\nStress: ");
                        for ( i = 1; i <= fullStressVector.giveSize(); i++ ) {
                            printf( " %e", fullStressVector.at(i) );
                        }

                        printf("\nKappa : ");
                        for ( i = 1; i <= strainSpaceHardeningVariables.giveSize(); i++ ) {
                            printf( " %e", strainSpaceHardeningVariables.at(i) );
                        }

                        printf("\n Restart %d\n", restart);
                    }

#endif

                    if ( restart ) {
                        IntArray newmap(nsurf);
                        actSurf = 0;
                        newmap.zero();
                        // buid new actsurf mask
                        for ( i = 1; i <= this->nsurf; i++ ) {
                            yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, strainSpaceHardeningVariables);
                            if ( ( ( activeConditionMap.at(i) ) && ( gamma.at(i) >= 0.0 ) ) || ( yieldValue > YIELD_TOL ) ) {
                                elastic = 0;
                                actSurf += 1;
                                newmap.at(i) = actSurf;
                                // if (actSurf == this->giveMaxNumberOfActiveYieldConds(gp)) break;
                            }
                        }

                        if ( this->getNewPopulation(activeConditionMap, newmap, this->giveMaxNumberOfActiveYieldConds(gp), nsurf) == 0 ) {
#ifdef DEBUG
                            if ( !debug ) {
                                this->clearPopulationSet();
                                debug = true;
                                goto huhu;
                            } else {
                                printf("ep:");
                                status->givePlasticStrainVector(plasticStrainVectorR);
                                plasticStrainVectorR.printYourself();
                                printf("kp:");
                                status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
                                strainSpaceHardeningVariables.printYourself();
                                elasticStrainVectorR = totalStrain;
                                elasticStrainVectorR.subtract(plasticStrainVectorR);
                                printf("sg:");
                                this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
                                fullStressVector.printYourself();

                                _error("Internal Consistency error: all combination of yield functions tryied, no consistent return");
                            }

#else
                            _error("Internal Consistency error: all combination of yield functions tryied, no consistent return");
#endif
                        }

                        for ( actSurf = 0, i = 1; i <= nsurf; i++ ) {
                            if ( activeConditionMap.at(i) ) {
                                actSurf++;
                            }
                        }

                        /*
                         * if (actSurf) {
                         * activeConditionMap=newmap;
                         * } else {
                         * restart = 0;
                         * // select any possible
                         * for (i=1; i<=this->nsurf; i++) {
                         *  if (initialConditionMap.at(i)) {
                         *    actSurf = 1;
                         *    activeConditionMap.zero();
                         *    activeConditionMap.at(i) = 1;
                         *    initialConditionMap.at(i) = 0;
                         *    restart = 1;
                         *    break;
                         *  }
                         * }
                         *
                         * if (!restart) {
                         *  _error ("Internal Consistency error");
                         * }
                         * }
                         */
                        status->givePlasticStrainVector(plasticStrainVectorR);
                        status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
                        elasticStrainVectorR = totalStrain;
                        elasticStrainVectorR.subtract(plasticStrainVectorR);
                        this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
                        gamma.zero();
                        nIterations = 0;
                        int restart = 0;
                        if ( restart ) {
                            goto huhu;
                        }

                        continue;
                    }

                    if ( yieldConsistency && ( residualVectorR.computeNorm() < RES_TOL ) ) {
                        return;
                    }
                }

                // compute  consistent tangent moduli
                this->computeAlgorithmicModuli(consistentModuli, gp, elasticModuliInverse,
                                               gamma, activeConditionMap,
                                               fullStressVector, strainSpaceHardeningVariables);


                /////////////////////////////////////////////////


                rhs.resize(actSurf);
                rhs.zero();
                gmat.resize(actSurf, actSurf);
                gmat.zero();
                lmat.resize(actSurf, strSize);
                lmat.zero();
                rmat.resize(strSize, actSurf);
                rmat.zero();
                if ( hasHardening ) {
                    this->computeReducedHardeningVarsSigmaGradient(ks, gp, activeConditionMap, fullStressVector,
                                                                   strainSpaceHardeningVariables, gamma);
                    this->computeReducedHardeningVarsLamGradient(kl, gp, actSurf, activeConditionMap,
                                                                 fullStressVector, strainSpaceHardeningVariables, gamma);
                }

                for ( i = 1; i <= nsurf; i++ ) {
                    if ( ( indx = activeConditionMap.at(i) ) ) {
                        if ( hasHardening ) {
                            helpVector.beTProductOf(ks, yieldGradKVec [ i - 1 ]);
                            helpVector.add(yieldGradSigVec [ i - 1 ]);
                            for ( j = 1; j <= strSize; j++ ) {
                                lmat.at(indx, j) = helpVector.at(j);
                            }
                        } else {
                            for ( j = 1; j <= strSize; j++ ) {
                                lmat.at(indx, j) = yieldGradSigVec [ i - 1 ].at(j);
                            }
                        }

                        if ( hasHardening ) {
                            this->computeReducedSKGradientMatrix(gradientMatrix,  i, gp, fullStressVector,
                                                                 strainSpaceHardeningVariables);
                            helpMtrx.beProductOf(gradientMatrix, kl);
                            helpMtrx.times( gamma.at(i) );
                            rmat.add(helpMtrx);
                        }

                        for ( j = 1; j <= strSize; j++ ) {
                            rmat.at(j, indx) += ( * loadGradSigVecPtr ) [ i - 1 ].at(j);
                        }

                        if ( hasHardening ) {
                            helpVector2.beTProductOf(kl, yieldGradKVec [ i - 1 ]);
                            for ( j = 1; j <= actSurf; j++ ) {
                                gmat.at(indx, j) = ( -1.0 ) * helpVector2.at(j);
                            }
                        }

                        yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, strainSpaceHardeningVariables);
                        rhs.at(indx) = yieldValue;
                        //helpVector2.beTProductOf (consistentModuli, helpVector);
                        //rhs.at(indx)=yieldValue-dotProduct(helpVector2, residualVectorR, residualVectorR.giveSize());
                    }
                }

                helpMtrx.beProductOf(lmat, consistentModuli);
                helpMtrx2.beProductOf(helpMtrx, rmat);
                gmat.add(helpMtrx2);
                helpVector2.beProductOf(consistentModuli, residualVectorR);
                helpVector.beProductOf(lmat, helpVector2);
                rhs.subtract(helpVector);


                // obtain increment to consistency parameter
                gmat.solveForRhs(rhs, dgammared);

                // assign gammas
                for ( i = 1; i <= this->nsurf; i++ ) {
                    if ( ( indx = activeConditionMap.at(i) ) ) {
                        dgamma.at(i) = dgammared.at(indx);
                    } else {
                        dgamma.at(i) = 0.0;
                    }
                }

                tempGamma = gamma;
                tempGamma.add(dgamma);


                if ( iterativeUpdateOfActiveConds ) {
                    //
                    // check active constraint
                    //

                    restart = 0;
                    for ( i = 1; i <= this->nsurf; i++ ) {
                        if ( tempGamma.at(i) < 0.0 ) {
                            restart = 1;
                        }
                    }

                    if ( restart ) {
                        // build up new activeConditionMap and restart the iterartion from begining
                        actSurf = 0;
                        activeConditionMap.zero();
                        for ( i = 1; i <= this->nsurf; i++ ) {
                            if ( tempGamma.at(i) > 0.0 ) {
                                actSurf += 1;
                                activeConditionMap.at(i) = actSurf;
                            }
                        }

                        status->givePlasticStrainVector(plasticStrainVectorR);
                        status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
                        elasticStrainVectorR = totalStrain;
                        elasticStrainVectorR.subtract(plasticStrainVectorR);
                        this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
                        gamma.zero();
                        nIterations = 0;
                        continue;
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            } while ( 1 ); // end restart loop

            // obtain incremental plastic strains and internal variables
            // residualVectorR vector and gradVec array are changed
            // but they are computed once again when new iteration starts

            helpVector2.resize(strSize);
            // Update state variables and consistency parameter
            for ( i = 1; i <= nsurf; i++ ) {
                if ( activeConditionMap.at(i) ) {
                    ( * loadGradSigVecPtr ) [ i - 1 ].times( dgamma.at(i) );
                    residualVectorR.add( ( * loadGradSigVecPtr ) [ i - 1 ] );

                    if ( hasHardening ) {
                        this->computeReducedSKGradientMatrix(gradientMatrix,  i, gp, fullStressVector,
                                                             strainSpaceHardeningVariables);
                        this->computeReducedHardeningVarsLamGradient(kl, gp, actSurf, activeConditionMap,
                                                                     fullStressVector, strainSpaceHardeningVariables, gamma);
                        helpMtrx.beProductOf(gradientMatrix, kl);
                        helpMtrx.times( gamma.at(i) );
                        helpVector2.beProductOf(helpMtrx, dgammared);

                        residualVectorR.add(helpVector2);
                    }
                }
            }

            helpVector.beProductOf(consistentModuli, residualVectorR);
            helpVector2.beProductOf(elasticModuliInverse, helpVector);

            gamma.add(dgamma);
            for ( i = 1; i <= strSize; i++ ) {
                plasticStrainVectorR.at(i) += helpVector2.at(i);
            }

            elasticStrainVectorR = totalStrain;
            elasticStrainVectorR.subtract(plasticStrainVectorR);
            // stress vector in full form due to computational convinience
            this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);

            if ( hasHardening ) {
                this->computeStrainHardeningVarsIncrement(dkappa, gp, fullStressVector, gamma, helpVector2, activeConditionMap);
                status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
                strainSpaceHardeningVariables.add(dkappa);
            }

            //this->computeStrainHardeningVarsIncrement(dkappa, fullStressVector, dgamma);
            //strainSpaceHardeningVariables.add(dkappa);
            //this-> computeTrialStressIncrement (fullStressVector, gp, elasticStrainVectorR, atTime);

            // increment iteration count
            nIterations++;

            if ( nIterations > PLASTIC_MATERIAL_MAX_ITERATIONS ) {
#ifdef DEBUG
                // issue warning
                fprintf( stderr, "\nclosestPointReturn (mat no. %d): reached max number of iterations\n", this->giveNumber() );
                fprintf(stderr, "activeSet: ");
                for ( i = 1; i <= nsurf; i++ ) {
                    fprintf( stderr, "%d ", activeConditionMap.at(i) );
                }

                fprintf(stderr, "\n");
                if ( !debug ) {
                    debug = true;
                    this->clearPopulationSet();
                    goto huhu;
                }

#endif


                // try last resort select single function from initial active set
                restart = 0;
                bool smartCandidate = false;
                int candidate = 0;
                IntArray newmap(nsurf);
                actSurf = 0;
                newmap.zero();
                // check for active conditions
                for ( i = 1; i <= this->nsurf; i++ ) {
                    if ( activeConditionMap.at(i) ) {
                        if ( gamma.at(i) < 0.0 ) {
                            smartCandidate = true;
                        } else {
                            candidate = i;
                        }
                    }
                }

                if ( smartCandidate && candidate && initialConditionMap.at(candidate) ) {
                    actSurf = 1;
                    IntArray newMap(nsurf);
                    newMap.at(candidate) = 1;

                    //initialConditionMap.at(candidate) = 0;
                    if ( this->getNewPopulation(activeConditionMap, newmap, this->giveMaxNumberOfActiveYieldConds(gp), nsurf) == 0 ) {
#ifdef DEBUG
                        if ( !debug ) {
                            this->clearPopulationSet();
                            debug = true;
                            goto huhu;
                        } else {
                            printf("ep:");
                            status->givePlasticStrainVector(plasticStrainVectorR);
                            plasticStrainVectorR.printYourself();
                            printf("kp:");
                            status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
                            strainSpaceHardeningVariables.printYourself();
                            elasticStrainVectorR = totalStrain;
                            elasticStrainVectorR.subtract(plasticStrainVectorR);
                            printf("sg:");
                            this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
                            fullStressVector.printYourself();

                            _error("Internal Consistency error: all combination of yield functions tryied, no consistent return");
                        }

#else
                        _error("Internal Consistency error: all combination of yield functions tryied, no consistent return");
#endif
                    }

                    for ( actSurf = 0, i = 1; i <= nsurf; i++ ) {
                        if ( activeConditionMap.at(i) ) {
                            actSurf++;
                        }
                    }

                    restart = 1;
                } else {
                    IntArray newMap(nsurf);
                    // select any possible
                    for ( i = 1; i <= this->nsurf; i++ ) {
                        if ( initialConditionMap.at(i) ) {
                            actSurf = 1;
                            newMap.at(i) = 1;
                            //initialConditionMap.at(i) = 0;
                            if ( this->getNewPopulation(activeConditionMap, newmap, this->giveMaxNumberOfActiveYieldConds(gp), nsurf) == 0 ) {
#ifdef DEBUG
                                if ( !debug ) {
                                    this->clearPopulationSet();
                                    debug = true;
                                    goto huhu;
                                } else {
                                    printf("ep:");
                                    status->givePlasticStrainVector(plasticStrainVectorR);
                                    plasticStrainVectorR.printYourself();
                                    printf("kp:");
                                    status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
                                    strainSpaceHardeningVariables.printYourself();
                                    elasticStrainVectorR = totalStrain;
                                    elasticStrainVectorR.subtract(plasticStrainVectorR);
                                    printf("sg:");
                                    this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
                                    fullStressVector.printYourself();

                                    _error("Internal Consistency error: all combination of yield functions tryied, no consistent return");
                                }

#else
                                _error("Internal Consistency error: all combination of yield functions tryied, no consistent return");
#endif
                            }

                            for ( actSurf = 0, i = 1; i <= nsurf; i++ ) {
                                if ( activeConditionMap.at(i) ) {
                                    actSurf++;
                                }
                            }

                            restart = 1;
                            break;
                        }
                    }
                }

                if ( restart ) {
                    //fprintf(stderr,"===>LAST RESORT<====");

                    status->givePlasticStrainVector(plasticStrainVectorR);
                    status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
                    elasticStrainVectorR = totalStrain;
                    elasticStrainVectorR.subtract(plasticStrainVectorR);
                    this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
                    gamma.zero();
                    nIterations = 0;
                    int restart = 0;
                    if ( restart ) {
                        goto huhu;
                    }

                    continue;
                } else {
                    _warning4( "GiveRealStressVector: local equlibrium not reached in %d iterations\nElement %d, gp %d, continuing",
                              PLASTIC_MATERIAL_MAX_ITERATIONS, gp->giveElement()->giveNumber(), gp->giveNumber() );
                    answer = fullStressVector;
                    // debug line
                    nIterations = 0;
                    goto huhu;
                    break;
                }
            }
        } while ( 1 );
    } else {
        answer = fullStressVector;
    }
}

void
MPlasticMaterial2 :: cuttingPlaneReturn(FloatArray &answer,
                                        IntArray &activeConditionMap,
                                        FloatArray &gamma,
                                        MatResponseForm form,
                                        GaussPoint *gp,
                                        const FloatArray &totalStrain,
                                        FloatArray &plasticStrainVectorR,
                                        FloatArray &strainSpaceHardeningVariables,
                                        TimeStep *atTime)
{
    FloatArray elasticStrainVectorR;
    FloatArray fullStressVector;
    FloatArray fSigmaGradientVectorR, gGradientVectorR, helpVector, helpVector2, rhs, dkappa;
    FloatArray di, dj;
    FloatMatrix elasticModuli, helpMtrx, helpMtrx2, gmat;
    FloatMatrix kl, ks, lmat, rmat;
    IntArray initialConditionMap;
    std :: vector< FloatArray >yieldGradSigVec(this->nsurf), loadGradSigVec(this->nsurf), * yieldGradSigVecPtr, * loadGradSigVecPtr;
    std :: vector< FloatArray >yieldGradKVec(this->nsurf), loadGradKVec(this->nsurf), * yieldGradKVecPtr, * loadGradKVecPtr;
    FloatArray dgamma(this->nsurf);
    double yieldValue;
    int nIterations = 0;
    int strSize, i, j, elastic, restart, actSurf, indx;
    bool yieldConsistency, init = true;
    bool hasHardening = this->hasHardening();

    MPlasticMaterial2Status *status = ( MPlasticMaterial2Status * ) this->giveStatus(gp);

    if ( this->plType == associatedPT ) {
        yieldGradSigVecPtr = loadGradSigVecPtr = & yieldGradSigVec;
        yieldGradKVecPtr = loadGradKVecPtr = & yieldGradKVec;
    } else {
        yieldGradSigVecPtr = & yieldGradSigVec;
        loadGradSigVecPtr  = & loadGradSigVec;
        yieldGradKVecPtr = & yieldGradKVec;
        loadGradKVecPtr = & loadGradKVec;
    }

huhu:
    status->givePlasticStrainVector(plasticStrainVectorR);
    status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);


    dgamma.resize(nsurf);
    dgamma.zero();
    gamma.resize(nsurf);
    gamma.zero();
    // compute elastic moduli
    this->computeReducedElasticModuli(elasticModuli, gp, atTime);
    strSize = elasticModuli.giveNumberOfRows();

    // initialize activeConditionMap
    elasticStrainVectorR = totalStrain;
    elasticStrainVectorR.subtract(plasticStrainVectorR);
    // stress vector in full form due to computational convinience
    this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);

    elastic = 1;
    actSurf = 0;
    activeConditionMap.zero();
    if ( init ) {
        initialConditionMap.resize(nsurf);
        initialConditionMap.zero();
    }

    for ( i = 1; i <= this->nsurf; i++ ) {
        yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, strainSpaceHardeningVariables);
        if ( yieldValue > YIELD_TOL ) {
            elastic = 0;
            if ( init ) {
                initialConditionMap.at(i) = 1;
            }

            if ( actSurf < this->giveMaxNumberOfActiveYieldConds(gp) ) {
                actSurf += 1;
                activeConditionMap.at(i) = actSurf;
            }
        }
    }

    init = false;

    if ( !elastic ) {
        do {            // local equilibrium loop
            do {              // loop for restart
                // compute gradients
                for ( i = 1; i <= nsurf; i++ ) {
                    this->computeReducedStressGradientVector(yieldGradSigVec [ i - 1 ], yieldFunction, i, gp, fullStressVector,
                                                             strainSpaceHardeningVariables);

                    if ( hasHardening ) {
                        this->computeKGradientVector(yieldGradKVec [ i - 1 ], yieldFunction, i, gp, fullStressVector,
                                                     strainSpaceHardeningVariables);
                    }
                }

                if ( this->plType == nonassociatedPT ) {
                    for ( i = 1; i <= nsurf; i++ ) {
                        this->computeReducedStressGradientVector(loadGradSigVec [ i - 1 ], loadFunction, i, gp, fullStressVector,
                                                                 strainSpaceHardeningVariables);
                        if ( hasHardening ) {
                            this->computeKGradientVector(loadGradKVec [ i - 1 ], loadFunction, i, gp, fullStressVector,
                                                         strainSpaceHardeningVariables);
                        }
                    }
                }

                rhs.resize(actSurf);
                rhs.zero();
                gmat.resize(actSurf, actSurf);
                gmat.zero();
                lmat.resize(actSurf, strSize);
                lmat.zero();
                rmat.resize(strSize, actSurf);
                rmat.zero();

                if ( hasHardening ) {
                    this->computeReducedHardeningVarsSigmaGradient(ks, gp, activeConditionMap, fullStressVector,
                                                                   strainSpaceHardeningVariables, gamma);
                    this->computeReducedHardeningVarsLamGradient(kl, gp, actSurf, activeConditionMap,
                                                                 fullStressVector, strainSpaceHardeningVariables, gamma);
                }

                for ( i = 1; i <= nsurf; i++ ) {
                    if ( ( indx = activeConditionMap.at(i) ) ) {
                        rhs.at(indx) = this->computeYieldValueAt(gp, i, fullStressVector, strainSpaceHardeningVariables);

                        if ( hasHardening ) {
                            helpVector.beTProductOf(ks, yieldGradKVec [ i - 1 ]);
                            helpVector.add(yieldGradSigVec [ i - 1 ]);
                            for ( j = 1; j <= strSize; j++ ) {
                                lmat.at(indx, j) = helpVector.at(j);
                            }
                        } else {
                            for ( j = 1; j <= strSize; j++ ) {
                                lmat.at(indx, j) = yieldGradSigVec [ i - 1 ].at(j);
                            }
                        }

                        for ( j = 1; j <= strSize; j++ ) {
                            rmat.at(j, indx) += ( * loadGradSigVecPtr ) [ i - 1 ].at(j);
                        }

                        if ( hasHardening ) {
                            helpVector2.beTProductOf(kl, yieldGradKVec [ i - 1 ]);
                            for ( j = 1; j <= actSurf; j++ ) {
                                gmat.at(indx, j) = ( -1.0 ) * helpVector2.at(j);
                            }
                        }
                    }
                }

                helpMtrx.beProductOf(lmat, elasticModuli);
                helpMtrx2.beProductOf(helpMtrx, rmat);
                gmat.add(helpMtrx2);

                // solve for plastic multiplier increments
                gmat.solveForRhs(rhs, helpVector);
                // assign gammas
                for ( i = 1; i <= this->nsurf; i++ ) {
                    if ( ( indx = activeConditionMap.at(i) ) ) {
                        dgamma.at(i) = helpVector.at(indx);
                    } else {
                        dgamma.at(i) = 0.0;
                    }
                }

                if ( iterativeUpdateOfActiveConds ) {
                    restart = 0;
                    // check for active conditions
                    for ( i = 1; i <= this->nsurf; i++ ) {
                        if ( ( gamma.at(i) + dgamma.at(i) ) < 0.0 ) {
                            gamma.at(i) = 0.0;
                            restart = 1;
                        }
                    }

                    if ( restart ) {
                        // build up new activeConditionMap and restart the iterartion from begining
                        actSurf = 0;
                        activeConditionMap.zero();
                        for ( i = 1; i <= this->nsurf; i++ ) {
                            if ( ( gamma.at(i) + dgamma.at(i) ) > 0.0 ) {
                                actSurf += 1;
                                activeConditionMap.at(i) = actSurf;
                            }
                        }

                        status->givePlasticStrainVector(plasticStrainVectorR);
                        status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
                        elasticStrainVectorR = totalStrain;
                        elasticStrainVectorR.subtract(plasticStrainVectorR);
                        this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
                        gamma.zero();
                        nIterations = 0;
                        continue;
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            } while ( 1 );           // end restart loop

            // compute plastic strain vector
            helpVector.resize(strSize);
            helpVector.zero();
            for ( i = 1; i <= this->nsurf; i++ ) {
                if ( activeConditionMap.at(i) ) {
                    ( * loadGradSigVecPtr ) [ i - 1 ].times( dgamma.at(i) );
                    for ( j = 1; j <= strSize; j++ ) {
                        helpVector.at(j) += ( * loadGradSigVecPtr ) [ i - 1 ].at(j);
                    }
                }
            }

            plasticStrainVectorR.add(helpVector);

            // update plastic multipliers
            gamma.add(dgamma);

            elasticStrainVectorR = totalStrain;
            elasticStrainVectorR.subtract(plasticStrainVectorR);
            // stress vector in full form due to computational convinience
            this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);

            if ( hasHardening ) {
                this->computeStrainHardeningVarsIncrement(dkappa, gp, fullStressVector, gamma, helpVector, activeConditionMap);
                status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
                strainSpaceHardeningVariables.add(dkappa);
            }

            // check convergence
            yieldConsistency = true;
            for ( i = 1; i <= nsurf; i++ ) {
                if ( activeConditionMap.at(i) ) {
                    yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, strainSpaceHardeningVariables);
                    if ( fabs(yieldValue) >= YIELD_TOL ) {
                        yieldConsistency = false;
                    }
                }
            }

            if ( yieldConsistency ) {
                restart = 0;
                // check for active conditions
                for ( i = 1; i <= this->nsurf; i++ ) {
                    yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, strainSpaceHardeningVariables);
                    if ( ( gamma.at(i) < 0.0 ) || ( ( activeConditionMap.at(i) == 0 ) && ( yieldValue > YIELD_TOL ) ) ) {
                        restart = 1;
                        //printf ("*");
                    }
                }

                if ( !restart ) {
                    break;
                } else {
                    IntArray newmap(nsurf);
                    actSurf = 0;
                    newmap.zero();
                    // buid new actsurf mask
                    for ( i = 1; i <= this->nsurf; i++ ) {
                        yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, strainSpaceHardeningVariables);
                        if ( ( ( activeConditionMap.at(i) ) && ( gamma.at(i) >= 0.0 ) ) || ( yieldValue > YIELD_TOL ) ) {
                            elastic = 0;
                            actSurf += 1;
                            newmap.at(i) = actSurf;
                            if ( actSurf == this->giveMaxNumberOfActiveYieldConds(gp) ) {
                                break;
                            }
                        }
                    }

                    activeConditionMap = newmap;

                    if ( actSurf ) {
                        activeConditionMap = newmap;
                    } else {
                        restart = 0;
                        // select any possible
                        for ( i = 1; i <= this->nsurf; i++ ) {
                            if ( initialConditionMap.at(i) ) {
                                actSurf = 1;
                                activeConditionMap.zero();
                                activeConditionMap.at(i) = 1;
                                initialConditionMap.at(i) = 0;
                                restart = 1;
                                break;
                            }
                        }
                    }




                    status->givePlasticStrainVector(plasticStrainVectorR);
                    status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
                    elasticStrainVectorR = totalStrain;
                    elasticStrainVectorR.subtract(plasticStrainVectorR);
                    this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
                    gamma.zero();
                    nIterations = 0;
                    int restart = 0;
                    if ( restart ) {
                        goto huhu;
                    }

                    continue;
                    //int restart = 0;
                    //if (restart) goto huhu;
                }
            }

            // increment iteration count
            nIterations++;

            if ( nIterations > PLASTIC_MATERIAL_MAX_ITERATIONS ) {
                // try last resort select single function from initial active set
                restart = 0;
                bool smartCandidate = false;
                int candidate = 0;
                IntArray newmap(nsurf);
                actSurf = 0;
                newmap.zero();
                // check for active conditions
                for ( i = 1; i <= this->nsurf; i++ ) {
                    if ( activeConditionMap.at(i) ) {
                        if ( gamma.at(i) < 0.0 ) {
                            smartCandidate = true;
                        } else {
                            candidate = i;
                        }
                    }
                }

                if ( smartCandidate && candidate && initialConditionMap.at(candidate) ) {
                    actSurf = 1;
                    activeConditionMap.zero();
                    activeConditionMap.at(candidate) = 1;
                    initialConditionMap.at(candidate) = 0;
                    restart = 1;
                } else {
                    // select any possible
                    for ( i = 1; i <= this->nsurf; i++ ) {
                        if ( initialConditionMap.at(i) ) {
                            actSurf = 1;
                            activeConditionMap.zero();
                            activeConditionMap.at(i) = 1;
                            initialConditionMap.at(i) = 0;
                            restart = 1;
                            break;
                        }
                    }
                }

                if ( restart ) {
                    //fprintf(stderr,"===>LAST RESORT<====");

                    status->givePlasticStrainVector(plasticStrainVectorR);
                    status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
                    elasticStrainVectorR = totalStrain;
                    elasticStrainVectorR.subtract(plasticStrainVectorR);
                    this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
                    gamma.zero();
                    nIterations = 0;
                    int restart = 0;
                    if ( restart ) {
                        goto huhu;
                    }

                    continue;
                } else {
                    _warning4( "GiveRealStressVector: local equlibrium not reached in %d iterations\nElement %d, gp %d, continuing",
                              PLASTIC_MATERIAL_MAX_ITERATIONS, gp->giveElement()->giveNumber(), gp->giveNumber() );
                    answer = fullStressVector;
                    // debug line
                    nIterations = 0;
                    goto huhu;
                    break;
                }
            }
        } while ( 1 );
    }

    //printf (" (%d iterations)", nIterations);
    answer = fullStressVector;
}



void
MPlasticMaterial2 :: computeReducedStressGradientVector(FloatArray &answer, functType ftype, int isurf,
                                                        GaussPoint *gp, const FloatArray &stressVector,
                                                        const FloatArray &strainSpaceHardeningVariables)
{
    FloatArray stressGradient;
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );

    this->computeStressGradientVector(stressGradient, ftype, isurf, gp, stressVector,
                                      strainSpaceHardeningVariables);

    crossSection->giveReducedCharacteristicVector(answer, gp, stressGradient);
}

void
MPlasticMaterial2 :: computeResidualVector(FloatArray &answer, GaussPoint *gp, const FloatArray &gamma,
                                           const IntArray &activeConditionMap, const FloatArray &plasticStrainVectorR,
                                           const FloatArray &strainSpaceHardeningVariables,
                                           std :: vector< FloatArray > &gradientVectorR)
{
    /* Computes Residual vector for closes point projection algorithm */

    FloatArray oldPlasticStrainVectorR;
    int i, j, size;
    MPlasticMaterial2Status *status = ( MPlasticMaterial2Status * ) this->giveStatus(gp);

    size = plasticStrainVectorR.giveSize();

    answer.resize(size);
    status->givePlasticStrainVector(oldPlasticStrainVectorR);

    for ( i = 1; i <= size; i++ ) {
        answer.at(i) = oldPlasticStrainVectorR.at(i) - plasticStrainVectorR.at(i);
    }

    for ( i = 1; i <= this->nsurf; i++ ) {
        if ( activeConditionMap.at(i) ) {
            for ( j = 1; j <= size; j++ ) {
                answer.at(j) += gamma.at(i) * gradientVectorR [ i - 1 ].at(j);
            }
        }
    }
}






void
MPlasticMaterial2 :: computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                                 const FloatArray &elasticStrainVectorR,
                                                 TimeStep *atTime)
{
    /* Computes the full trial elastic stress vector */

    FloatMatrix de;
    FloatArray reducedAnswer;
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );

    this->computeReducedElasticModuli(de, gp, atTime);
    /*
     *    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(de, ReducedForm, TangentStiffness,
     *                                                                                                                                                                                                                                                    gp, atTime) ;
     */
    reducedAnswer.beProductOf(de, elasticStrainVectorR);
    crossSection->giveFullCharacteristicVector(answer, gp, reducedAnswer);
}


void
MPlasticMaterial2 :: computeAlgorithmicModuli(FloatMatrix &answer,
                                              GaussPoint *gp,
                                              const FloatMatrix &elasticModuliInverse,
                                              const FloatArray &gamma,
                                              const IntArray &activeConditionMap,
                                              const FloatArray &fullStressVector,
                                              const FloatArray &strainSpaceHardeningVariables)
{
    /* returns consistent moduli in reduced form.
     * stressVector in full form
     */
    FloatMatrix gradientMatrix, ks;
    FloatMatrix help, helpInverse;
    int i, j, size;
    bool hasHardening = this->hasHardening();

    size = elasticModuliInverse.giveNumberOfRows();

    // assemble consistent moduli
    helpInverse.resize(size, size);
    helpInverse.zero();
    for ( i = 1; i <= nsurf; i++ ) {
        if ( activeConditionMap.at(i) ) {
            // ask gradientMatrix
            this->computeReducedSSGradientMatrix(gradientMatrix,  i, gp, fullStressVector,
                                                 strainSpaceHardeningVariables);

            gradientMatrix.times( gamma.at(i) );
            helpInverse.add(gradientMatrix);

            if ( hasHardening ) {
                this->computeReducedSKGradientMatrix(gradientMatrix,  i, gp, fullStressVector,
                                                     strainSpaceHardeningVariables);

                gradientMatrix.times( gamma.at(i) );
                this->computeReducedHardeningVarsSigmaGradient(ks, gp, activeConditionMap, fullStressVector, strainSpaceHardeningVariables, gamma);
                help.beProductOf(gradientMatrix, ks);
                helpInverse.add(help);
            }
        }
    }

    for ( i = 1; i <= size; i++ ) {
        for ( j = 1; j <= size; j++ ) {
            helpInverse.at(i, j) += elasticModuliInverse.at(i, j);
        }
    }

    // prevent too big numbers
    for ( i = 1; i <= size; i++ ) {
        if ( fabs( helpInverse.at(i, i) ) > 1.e16 ) {
            helpInverse.at(i, i) = sgn( helpInverse.at(i, i) ) * 1.e16;
        }
    }

    answer.beInverseOf(helpInverse);
}

// ----------------------------------------------------------------------------//


void
MPlasticMaterial2 :: giveConsistentStiffnessMatrix(FloatMatrix &answer,
                                                   MatResponseForm form,
                                                   MatResponseMode mode,
                                                   GaussPoint *gp,
                                                   TimeStep *atTime)
{
    //
    // returns receiver material matrix for given reached state
    // described by temp variables.
    //

    FloatMatrix consistentModuli, elasticModuli, umat, vmat;
    FloatMatrix consistentModuliInverse, elasticModuliInverse, kl, ks;
    FloatMatrix gradientMatrix, gmat, gmatInv, gradMat, helpMtrx, helpMtrx2, answerR;
    FloatArray gradientVector, stressVector, fullStressVector;
    FloatArray strainSpaceHardeningVariables, helpVector;
    std :: vector< FloatArray >yieldGradSigVec(this->nsurf), loadGradSigVec(this->nsurf), * yieldGradSigVecPtr, * loadGradSigVecPtr;
    std :: vector< FloatArray >yieldGradKVec(this->nsurf), loadGradKVec(this->nsurf), * yieldGradKVecPtr, * loadGradKVecPtr;
    FloatArray helpVector2;

    IntArray activeConditionMap, mask;
    FloatArray gamma;
    int strSize, i, j, indx, actSurf = 0;
    bool hasHardening = this->hasHardening();

    MPlasticMaterial2Status *status = ( MPlasticMaterial2Status * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );

    if ( this->plType == associatedPT ) {
        yieldGradSigVecPtr = loadGradSigVecPtr = & yieldGradSigVec;
        yieldGradKVecPtr = loadGradKVecPtr = & yieldGradKVec;
    } else {
        yieldGradSigVecPtr = & yieldGradSigVec;
        loadGradSigVecPtr  = & loadGradSigVec;
        yieldGradKVecPtr = & yieldGradKVec;
        loadGradKVecPtr = & loadGradKVec;
    }

    // ask for plastic consistency parameter
    status->giveTempGamma(gamma);
    status->giveTempActiveConditionMap(activeConditionMap);
    //
    // check for elastic cases
    //
    if ( ( status->giveTempStateFlag() == MPlasticMaterial2Status :: PM_Elastic ) ||
        ( status->giveTempStateFlag() == MPlasticMaterial2Status :: PM_Unloading ) ) {
        this->giveCharacteristicMatrix(answer, form, ElasticStiffness, gp, atTime);
        return;
    }

    //
    // plastic case
    //
    // determine number of active surfaces
    for ( i = 1; i <= nsurf; i++ ) {
        if ( activeConditionMap.at(i) ) {
            actSurf++;
        }
    }

    // compute elastic moduli and its inverse
    this->computeReducedElasticModuli(elasticModuli, gp, atTime);
    elasticModuliInverse.beInverseOf(elasticModuli);
    strSize = elasticModuliInverse.giveNumberOfRows();

    stressVector = status->giveTempStressVector();
    crossSection->giveFullCharacteristicVector(fullStressVector, gp, stressVector);
    status->giveTempStrainSpaceHardeningVarsVector(strainSpaceHardeningVariables);

    //
    // compute consistent moduli
    //
    this->computeAlgorithmicModuli(consistentModuli, gp, elasticModuliInverse, gamma,
                                   activeConditionMap, fullStressVector, strainSpaceHardeningVariables);

    //computee gmatInv
    for ( i = 1; i <= nsurf; i++ ) {
        this->computeReducedStressGradientVector(yieldGradSigVec [ i - 1 ], yieldFunction, i, gp, fullStressVector,
                                                 strainSpaceHardeningVariables);
        if ( hasHardening ) {
            this->computeKGradientVector(yieldGradKVec [ i - 1 ], yieldFunction, i, gp, fullStressVector,
                                         strainSpaceHardeningVariables);
        }
    }

    if ( this->plType == nonassociatedPT ) {
        for ( i = 1; i <= nsurf; i++ ) {
            this->computeReducedStressGradientVector(loadGradSigVec [ i - 1 ], loadFunction, i, gp, fullStressVector,
                                                     strainSpaceHardeningVariables);
            if ( hasHardening ) {
                this->computeKGradientVector(loadGradKVec [ i - 1 ], loadFunction, i, gp, fullStressVector,
                                             strainSpaceHardeningVariables);
            }
        }
    }


    umat.resize(strSize, actSurf);
    umat.zero();
    vmat.resize(actSurf, strSize);
    vmat.zero();
    gmat.resize(actSurf, actSurf);
    gmat.zero();

    if ( hasHardening ) {
        this->computeReducedHardeningVarsSigmaGradient(ks, gp, activeConditionMap, fullStressVector,
                                                       strainSpaceHardeningVariables, gamma);
        this->computeReducedHardeningVarsLamGradient(kl, gp, actSurf, activeConditionMap,
                                                     fullStressVector, strainSpaceHardeningVariables, gamma);
    }

    for ( i = 1; i <= nsurf; i++ ) {
        if ( ( indx = activeConditionMap.at(i) ) ) {
            if ( hasHardening ) {
                helpVector.beTProductOf(ks, yieldGradKVec [ i - 1 ]);
                helpVector.add(yieldGradSigVec [ i - 1 ]);
                for ( j = 1; j <= strSize; j++ ) {
                    vmat.at(indx, j) = helpVector.at(j);
                }
            } else {
                for ( j = 1; j <= strSize; j++ ) {
                    vmat.at(indx, j) = yieldGradSigVec [ i - 1 ].at(j);
                }
            }

            if ( hasHardening ) {
                this->computeReducedSKGradientMatrix(gradientMatrix,  i, gp, fullStressVector,
                                                     strainSpaceHardeningVariables);
                helpMtrx.beProductOf(gradientMatrix, kl);
                helpMtrx.times( gamma.at(i) );
                umat.add(helpMtrx);
            }

            for ( j = 1; j <= strSize; j++ ) {
                umat.at(j, indx) += ( ( * loadGradSigVecPtr ) [ i - 1 ] ).at(j);
            }


            if ( hasHardening ) {
                helpVector2.beTProductOf(kl, yieldGradKVec [ i - 1 ]);
                for ( j = 1; j <= actSurf; j++ ) {
                    gmat.at(indx, j) = ( -1.0 ) * helpVector2.at(j);
                }
            }
        }
    }

    helpMtrx.beProductOf(vmat, consistentModuli); // V S
    helpMtrx2.beProductOf(helpMtrx, umat);
    gmat.add(helpMtrx2);
    /////////////////////////////
    gmatInv.beInverseOf(gmat);


    helpMtrx2.beProductOf(gmatInv, helpMtrx);
    helpMtrx.beProductOf(consistentModuli, umat); // S U
    answer.beProductOf(helpMtrx, helpMtrx2); // SUGVS
    answer.negated();
    answer.add(consistentModuli);

    if ( form == ReducedForm ) {
        return;
    } else {
        this->giveStressStrainMask( mask, ReducedForm, gp->giveMaterialMode() );
        answerR = answer;
        answer.beSubMatrixOfSizeOf(answerR, mask, 6);
    }
}

void
MPlasticMaterial2 :: giveElastoPlasticStiffnessMatrix(FloatMatrix &answer,
                                                      MatResponseForm form,
                                                      MatResponseMode mode,
                                                      GaussPoint *gp,
                                                      TimeStep *atTime)
{
    //
    // returns receiver material matrix for given reached state
    // described by temp variables.
    //

    FloatMatrix dmat, elasticModuli, umat, vmat;
    FloatMatrix gmat, gmatInv, helpMtrx, helpMtrx2, kl, ks;
    FloatArray gradientVector, stressVector, fullStressVector;
    FloatArray strainSpaceHardeningVariables, helpVector, helpVector2;
    std :: vector< FloatArray >yieldGradSigVec(this->nsurf), loadGradSigVec(this->nsurf), * yieldGradSigVecPtr, * loadGradSigVecPtr;
    std :: vector< FloatArray >yieldGradKVec(this->nsurf), loadGradKVec(this->nsurf), * yieldGradKVecPtr, * loadGradKVecPtr;
    FloatArray helpVec;

    IntArray activeConditionMap, mask;
    FloatArray gamma;
    int strSize, i, j, indx, actSurf = 0;
    bool hasHardening = this->hasHardening();

    MPlasticMaterial2Status *status = ( MPlasticMaterial2Status * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );

    if ( this->plType == associatedPT ) {
        yieldGradSigVecPtr = loadGradSigVecPtr = & yieldGradSigVec;
        yieldGradKVecPtr = loadGradKVecPtr = & yieldGradKVec;
    } else {
        yieldGradSigVecPtr = & yieldGradSigVec;
        loadGradSigVecPtr  = & loadGradSigVec;
        yieldGradKVecPtr = & yieldGradKVec;
        loadGradKVecPtr = & loadGradKVec;
    }

    // ask for plastic consistency parameter
    status->giveTempGamma(gamma);
    status->giveTempActiveConditionMap(activeConditionMap);
    //
    // check for elastic cases
    //
    if ( ( status->giveTempStateFlag() == MPlasticMaterial2Status :: PM_Elastic ) ||
        ( status->giveTempStateFlag() == MPlasticMaterial2Status :: PM_Unloading ) ) {
        this->giveCharacteristicMatrix(answer, form, ElasticStiffness, gp, atTime);
        return;
    }

    //
    // plastic case
    //
    // determine number of active surfaces
    for ( i = 1; i <= nsurf; i++ ) {
        if ( activeConditionMap.at(i) ) {
            actSurf++;
        }
    }

    // compute elastic moduli and its inverse
    this->computeReducedElasticModuli(elasticModuli, gp, atTime);


    stressVector = status->giveStressVector();
    crossSection->giveFullCharacteristicVector(fullStressVector, gp, stressVector);
    status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
    strSize = elasticModuli.giveNumberOfRows();

    //computee gmatInv
    for ( i = 1; i <= nsurf; i++ ) {
        this->computeReducedStressGradientVector(yieldGradSigVec [ i - 1 ], yieldFunction, i, gp, fullStressVector,
                                                 strainSpaceHardeningVariables);
        if ( hasHardening ) {
            this->computeKGradientVector(yieldGradKVec [ i - 1 ], yieldFunction, i, gp, fullStressVector,
                                         strainSpaceHardeningVariables);
        }
    }

    if ( this->plType == nonassociatedPT ) {
        for ( i = 1; i <= nsurf; i++ ) {
            this->computeReducedStressGradientVector(loadGradSigVec [ i - 1 ], loadFunction, i, gp, fullStressVector,
                                                     strainSpaceHardeningVariables);
            if ( hasHardening ) {
                this->computeKGradientVector(loadGradKVec [ i - 1 ], loadFunction, i, gp, fullStressVector,
                                             strainSpaceHardeningVariables);
            }
        }
    }

    umat.resize(strSize, actSurf);
    vmat.resize(actSurf, strSize);
    gmat.resize(actSurf, actSurf);

    if ( hasHardening ) {
        this->computeReducedHardeningVarsSigmaGradient(ks, gp, activeConditionMap, fullStressVector,
                                                       strainSpaceHardeningVariables, gamma);
        this->computeReducedHardeningVarsLamGradient(kl, gp, actSurf, activeConditionMap,
                                                     fullStressVector, strainSpaceHardeningVariables, gamma);
    }

    for ( i = 1; i <= nsurf; i++ ) {
        if ( ( indx = activeConditionMap.at(i) ) ) {
            if ( hasHardening ) {
                helpVector.beTProductOf(ks, yieldGradKVec [ i - 1 ]);
                helpVector.add(yieldGradSigVec [ i - 1 ]);
                for ( j = 1; j <= strSize; j++ ) {
                    vmat.at(indx, j) = helpVector.at(j);
                }
            } else {
                for ( j = 1; j <= strSize; j++ ) {
                    vmat.at(indx, j) = yieldGradSigVec [ i - 1 ].at(j);
                }
            }

            for ( j = 1; j <= strSize; j++ ) {
                umat.at(j, indx) += ( * loadGradSigVecPtr ) [ i - 1 ].at(j);
            }


            if ( hasHardening ) {
                helpVector2.beTProductOf(kl, yieldGradKVec [ i - 1 ]);
                for ( j = 1; j <= actSurf; j++ ) {
                    gmat.at(indx, j) = ( -1.0 ) * helpVector2.at(j);
                }
            }
        }
    }

    helpMtrx.beProductOf(vmat, elasticModuli); // V S
    helpMtrx2.beProductOf(helpMtrx, umat);
    gmat.add(helpMtrx2);
    /////////////////////////////
    gmatInv.beInverseOf(gmat);

    helpMtrx.beProductOf(gmatInv, vmat);
    helpMtrx2.beProductOf(umat, helpMtrx);
    answer.beProductOf(elasticModuli, helpMtrx2);
    answer.negated();
    answer.add(elasticModuli);

    if ( form == ReducedForm ) {
        return;
    } else {
        this->giveStressStrainMask( mask, ReducedForm, gp->giveMaterialMode() );
        helpMtrx = answer;
        answer.beSubMatrixOfSizeOf(helpMtrx, mask, 6);
    }
}

/*
 * void
 * MPlasticMaterial2 :: computeDiagModuli(FloatMatrix& answer,
 *                                    GaussPoint *gp, FloatMatrix &elasticModuliInverse,
 *                                    FloatMatrix &hardeningModuliInverse)
 * {
 * //
 * // assembles diagonal moduli from elasticModuliInverse and hardeningModuliInverse
 * //
 * int size1, size2, i,j;
 *
 * size1 = elasticModuliInverse.giveNumberOfRows();
 * if (hardeningModuliInverse.giveNumberOfRows()) size2=size1+hardeningModuliInverse.giveNumberOfRows();
 * else size2 = size1;
 *
 * answer.resize (size2, size2);
 * answer.zero();
 *
 * for (i=1; i<=size1; i++)
 * for (j=1; j<= size1; j++) answer.at(i,j) = elasticModuliInverse.at(i,j);
 *
 * for (i=size1+1; i<= size2; i++)
 * for (j=size1+1; j<= size2; j++) answer.at(i,j) = hardeningModuliInverse.at(i-size1, j-size1);
 * }
 */

void
MPlasticMaterial2 :: computeReducedElasticModuli(FloatMatrix &answer,
                                                 GaussPoint *gp,
                                                 TimeStep *atTime)
{  /* Returns elastic moduli in reduced stress-strain space*/
    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, ReducedForm,
                                                                ElasticStiffness,
                                                                gp, atTime);
}


// overloaded from structural material

void
MPlasticMaterial2 :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form,
                                                   MatResponseMode mode,
                                                   GaussPoint *gp,
                                                   TimeStep *atTime)
//
//
//
// computes full 3d constitutive matrix for case of 3d stress-strain state.
// it returns elasto-plastic stiffness material matrix.
// if strainIncrement == NULL a loading is assumed
// for detailed description see (W.F.Chen: Plasticity in Reinforced Concrete, McGraw-Hill, 1982,
// chapter 6.)
//
// if derived material would like to implement failure behaviour
// it must redefine basic Give3dMaterialStiffnessMatrix function
// in order to take possible failure (tension cracking) into account
//
//
//
{
    MaterialMode originalMode = gp->giveMaterialMode();
    if ( originalMode != _3dMat ) {
        _error("give3dMaterialStiffnessMatrix : Different stressStrain mode encountered");
    }

    // we can force 3d response, and we obtain correct 3d tangent matrix,
    // but in fact, stress integration algorithm will not work
    // because in stress integration algorithm we are unable to recognize
    // which reduction from 3d case should be performed to obtain correct result.
    // so for new stressStrain state, instead of programming 3d reduction,
    // you should enhance imposeConstraints functions for ne state, and
    // then programming simple inteface function for you stressstrain state
    // calling GiveMaterailStiffenssMatrix, which imposes constrains correctly.
    if ( mode == TangentStiffness ) {
        if ( rmType == mpm_ClosestPoint ) {
            this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
        } else {
            this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
        }
    } else {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    }
}


void
MPlasticMaterial2 :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form,
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
    if ( mode == TangentStiffness ) {
        if ( rmType == mpm_ClosestPoint ) {
            this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
        } else {
            this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
        }
    } else {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    }
}


void
MPlasticMaterial2 :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                              MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *atTime)

//
// return receiver's 2dPlaneStrainMtrx constructed from
// general 3dMatrialStiffnessMatrix
// (2dPlaneStrain ==> eps_z = gamma_xz = gamma_yz = 0.)
//
{
    if ( mode == TangentStiffness ) {
        if ( rmType == mpm_ClosestPoint ) {
            this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
        } else {
            this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
        }
    } else {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    }
}


void
MPlasticMaterial2 :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime)

//
// returns receiver's 1dMaterialStiffnessMAtrix
// (1d case ==> sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.)
{
    if ( mode == TangentStiffness ) {
        if ( rmType == mpm_ClosestPoint ) {
            this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
        } else {
            this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
        }
    } else {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    }
}



void
MPlasticMaterial2 :: give2dBeamLayerStiffMtrx(FloatMatrix &answer, MatResponseForm form,
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
    if ( mode == TangentStiffness ) {
        if ( rmType == mpm_ClosestPoint ) {
            this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
        } else {
            this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
        }
    } else {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    }
}



void
MPlasticMaterial2 :: give2dPlateLayerStiffMtrx(FloatMatrix &answer,
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
    if ( mode == TangentStiffness ) {
        if ( rmType == mpm_ClosestPoint ) {
            this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
        } else {
            this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
        }
    } else {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    }
}

void
MPlasticMaterial2 :: give1dFiberStiffMtrx(FloatMatrix &answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *atTime)
//
// returns receiver's 1dFiber
// (1dFiber ==> sigma_y = sigma_z = tau_yz = 0.)
//
// standard method from Material Class overloaded, because no inversion is needed.
// the reduction from 3d case will not work
// this implementation should be faster.
{
    if ( mode == TangentStiffness ) {
        if ( rmType == mpm_ClosestPoint ) {
            this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
        } else {
            this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
        }
    } else {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    }
}


void
MPlasticMaterial2 :: give3dShellLayerStiffMtrx(FloatMatrix &answer, MatResponseForm form,
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
    if ( mode == TangentStiffness ) {
        if ( rmType == mpm_ClosestPoint ) {
            this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
        } else {
            this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
        }
    } else {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    }
}


int
MPlasticMaterial2 :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    MPlasticMaterial2Status *status = ( MPlasticMaterial2Status * ) this->giveStatus(aGaussPoint);
    if ( type == IST_PlasticStrainTensor ) {
        status->givePlasticStrainVector(answer);
        return 1;
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        int indx;
        FloatArray st(6), s;

        status->givePlasticStrainVector(s);

        for ( int i = 1; i <= s.giveSize(); i++ ) {
            indx = this->giveStressStrainComponentIndOf(ReducedForm, aGaussPoint->giveMaterialMode(), i);
            if ( indx ) {
                st.at(indx) = s.at(i);
            }
        }

        this->computePrincipalValues(answer, st, principal_strain);
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}




InternalStateValueType
MPlasticMaterial2 :: giveIPValueType(InternalStateType type)
{
    // strains components packed in engineering notation
    if ( type == IST_PlasticStrainTensor ) {
        return ISVT_TENSOR_S3E;
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        return ISVT_VECTOR;
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }
}


int
MPlasticMaterial2 :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( type == IST_PlasticStrainTensor ) {
        this->giveStressStrainMask(answer, FullForm, mmode);
        return 1;
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        answer.resize(6);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        return 1;
    } else {
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}


int
MPlasticMaterial2 :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( type == IST_PlasticStrainTensor ) {
        return this->giveSizeOfReducedStressStrainVector( aGaussPoint->giveMaterialMode() );
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        return 3;
    } else {
        return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}

long
MPlasticMaterial2 :: getPopulationSignature(IntArray &mask)
{
    long val = 0;
    int size = mask.giveSize();
    if ( size > ( int ) ( sizeof( long ) * 8 ) ) {
        printf("MPlasticMaterial2: too many yield conditions");
    }

    for ( int i = 1; i <= mask.giveSize(); i++ ) {
        if ( mask.at(i) ) {
            val |= ( 1L << ( i - 1 ) );
        }
    }

    return val;
}

int
MPlasticMaterial2 :: testPopulation(long pop)
{
    if ( populationSet.find(pop) != populationSet.end() ) {
        return 0;
    } else {
        return 1;
    }
}

void
MPlasticMaterial2 :: clearPopulationSet()
{
    populationSet.clear();
}

void
MPlasticMaterial2 :: addNewPopulation(IntArray &mask)
{
    long val = getPopulationSignature(mask);
    populationSet.insert(val);
}


int
MPlasticMaterial2 :: getNewPopulation(IntArray &result, IntArray &candidateMask, int degree, int size)
{
    long val = 0;
    std :: set< long > :: iterator it;
    int i, candDegree = 0, candSize = candidateMask.giveSize();
    for ( i = 1; i <= candSize; i++ ) {
        if ( candidateMask.at(i) ) {
            candDegree++;
        }
    }

    result.resize(size);

    // first try candidate if compatible
    if ( candDegree <= degree ) {
        // generate candite popul value
        val = getPopulationSignature(candidateMask);

        if ( testPopulation(val) ) {
            // candite is o.k.
            populationSet.insert(val);
            result = candidateMask;
            return 1;
        }
    } else {
        // try to generate suitable candidate from candidateMask if possible
        // first construct array of candidates
        int j = 1;
        IntArray candidateArray(candDegree);
        IntArray ind(degree);
        for ( i = 1; i <= candSize; i++ ) {
            if ( candidateMask.at(i) ) {
                candidateArray.at(j++) = i;
            }
        }

        // first try first degree candidates
        int activeIndx = degree;
        for ( i = 1; i <= degree; i++ ) {
            ind.at(i) = i;
        }


        do {
            result.zero();
            int ii, count = 1;
            for ( i = 1; i <= degree; i++ ) {
                ii = candidateArray.at( ind.at(i) );
                if ( !result.at(ii) ) {
                    result.at(ii) = count++;
                }
            }

            // order result array
            count = 1;
            for ( i = 1; i <= size; i++ ) {
                if ( result.at(i) ) {
                    result.at(i) = count++;
                }
            }

            val = getPopulationSignature(result);
            if ( testPopulation(val) ) {
                // candite is o.k.
                populationSet.insert(val);
                return 1;
            }

            bool cont = false;
            activeIndx = degree;
            do {
                if ( ind.at(activeIndx) < candDegree ) {
                    ind.at(activeIndx)++;
                    cont = false;
                } else {
                    ind.at(activeIndx) = 1;
                    activeIndx--;
                    if ( activeIndx > 0 ) {
                        cont = true;
                    } else {
                        cont = false;
                    }
                }
            } while ( cont );
        } while ( activeIndx >= 1 );
    }

    // generate any suitable population
    IntArray ind(degree);
    // first try first degree candidates
    int activeIndx = degree;
    for ( i = 1; i <= degree; i++ ) {
        ind.at(i) = 1;
    }


    do {
        result.zero();
        int ii, count = 1;
        for ( i = 1; i <= degree; i++ ) {
            ii = ind.at(i);
            if ( !result.at(ii) ) {
                result.at(ii) = count++;
            }
        }

        // order result array
        count = 1;
        for ( i = 1; i <= size; i++ ) {
            if ( result.at(i) ) {
                result.at(i) = count++;
            }
        }

        val = getPopulationSignature(result);
        if ( testPopulation(val) ) {
            // candite is o.k.
            populationSet.insert(val);
            return 1;
        }

        bool cont = false;
        activeIndx = degree;
        do {
            if ( ind.at(activeIndx) < size ) {
                ind.at(activeIndx)++;
                cont = false;
            } else {
                ind.at(activeIndx) = 1;
                activeIndx--;
                if ( activeIndx > 0 ) {
                    cont = true;
                } else {
                    cont = false;
                }
            }
        }  while ( cont );
    } while ( activeIndx >= 1 );

    // there is no population left -> return 0;
    return 0;
}






#ifdef __OOFEG
#endif






//
// SmearedCrackingMaterialStatus Class
//

MPlasticMaterial2Status :: MPlasticMaterial2Status(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), plasticStrainVector(), tempPlasticStrainVector(),
    strainSpaceHardeningVarsVector(), tempStrainSpaceHardeningVarsVector()
{
    state_flag = temp_state_flag = MPlasticMaterial2Status :: PM_Elastic;
    gamma = tempGamma = 0.;
}

MPlasticMaterial2Status :: ~MPlasticMaterial2Status()
{ }

void
MPlasticMaterial2Status :: printOutputAt(FILE *file, TimeStep *tStep)
{
    int i, n;

    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( ( state_flag == MPlasticMaterial2Status :: PM_Yielding ) || ( state_flag == MPlasticMaterial2Status :: PM_Unloading ) ) {
        if ( state_flag == MPlasticMaterial2Status :: PM_Yielding ) {
            fprintf(file, " Yielding, ");
        } else {
            fprintf(file, " Unloading, ");
        }

        n = plasticStrainVector.giveSize();
        fprintf(file, " plastic strains ");
        for ( i = 1; i <= n; i++ ) {
            fprintf( file, " % .4e", plasticStrainVector.at(i) );
        }

        if ( strainSpaceHardeningVarsVector.giveSize() ) {
            n = strainSpaceHardeningVarsVector.giveSize();
            fprintf(file, ", strain space hardening vars ");
            for ( i = 1; i <= n; i++ ) {
                fprintf( file, " % .4e", strainSpaceHardeningVarsVector.at(i) );
            }
        }
    }

    fprintf(file, "}\n");
}


void MPlasticMaterial2Status :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
//
{
    StructuralMaterialStatus :: initTempStatus();

    if ( plasticStrainVector.giveSize() == 0 ) {
        plasticStrainVector.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->
                                   giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
        plasticStrainVector.zero();
    }

    tempPlasticStrainVector = plasticStrainVector;

    if ( strainSpaceHardeningVarsVector.giveSize() == 0 ) {
        strainSpaceHardeningVarsVector.resize( ( ( MPlasticMaterial2 * ) gp->giveMaterial() )->
                                              giveSizeOfReducedHardeningVarsVector(gp) );
        strainSpaceHardeningVarsVector.zero();
    }

    tempStrainSpaceHardeningVarsVector = strainSpaceHardeningVarsVector;

    temp_state_flag = state_flag;

    tempGamma = gamma;
    tempActiveConditionMap = activeConditionMap;
}



void
MPlasticMaterial2Status :: updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    StructuralMaterialStatus :: updateYourself(atTime);

    plasticStrainVector = tempPlasticStrainVector;
    strainSpaceHardeningVarsVector = tempStrainSpaceHardeningVarsVector;

    state_flag = temp_state_flag;
    gamma = tempGamma;
    activeConditionMap = tempActiveConditionMap;
}




contextIOResultType
MPlasticMaterial2Status :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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
    if ( ( iores = plasticStrainVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = strainSpaceHardeningVarsVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->write(& state_flag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = gamma.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = activeConditionMap.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    return CIO_OK;
}



contextIOResultType
MPlasticMaterial2Status :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = plasticStrainVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = strainSpaceHardeningVarsVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->read(& state_flag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = gamma.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = activeConditionMap.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    return CIO_OK; // return succes
}
} // end namespace oofem
