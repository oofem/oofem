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

#include "mplasticmaterial.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "structuralcrosssection.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
#define YIELD_TOL 1.e-4
#define RES_TOL   1.e-4
#define PLASTIC_MATERIAL_MAX_ITERATIONS 90

MPlasticMaterial :: MPlasticMaterial(int n, Domain *d)  : StructuralMaterial(n, d)
    //
    // constructor
    //
{
    linearElasticMaterial = NULL;
    rmType = mpm_ClosestPoint;
    plType = associatedPT;
}


MPlasticMaterial :: ~MPlasticMaterial()
//
// destructor
//
{
    delete linearElasticMaterial;
}


int
MPlasticMaterial :: hasMaterialModeCapability(MaterialMode mode)
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
MPlasticMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates new  material status  corresponding to this class
 */
{
    MPlasticMaterialStatus *status;

    status = new MPlasticMaterialStatus(1, this->giveDomain(), gp);
    return status;
}


void
MPlasticMaterial :: giveRealStressVector(FloatArray &answer,
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

    MPlasticMaterialStatus *status = ( MPlasticMaterialStatus * ) this->giveStatus(gp);
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
        newState = MPlasticMaterialStatus :: PM_Yielding;         // test if plastic loading occur
    }
    // no plastic loading - check for unloading
    else if ( ( state == MPlasticMaterialStatus :: PM_Yielding ) || ( state == MPlasticMaterialStatus :: PM_Unloading ) ) {
        newState = MPlasticMaterialStatus :: PM_Unloading;
    } else {
        newState = MPlasticMaterialStatus :: PM_Elastic;
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
MPlasticMaterial :: closestPointReturn(FloatArray &answer,
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

    FloatArray helpVector, helpVector2;
    FloatArray dgamma, tempGamma;
    FloatMatrix elasticModuli, hardeningModuli, consistentModuli;
    FloatMatrix elasticModuliInverse, hardeningModuliInverse;
    FloatMatrix helpMtrx, helpMtrx2;
    FloatMatrix gmat;
    std :: vector< FloatArray >yieldGradVec(this->nsurf), loadGradVec(this->nsurf), * yieldGradVecPtr, * loadGradVecPtr;
    FloatArray rhs;
    double yieldValue;
    int nIterations = 0;
    int i, j, strSize, totSize;
    int elastic, restart, actSurf, indx;
    bool yieldConsistency;


    if ( this->plType == associatedPT ) {
        yieldGradVecPtr = loadGradVecPtr = & yieldGradVec;
    } else {
        yieldGradVecPtr = & yieldGradVec;
        loadGradVecPtr  = & loadGradVec;
    }

    MPlasticMaterialStatus *status = ( MPlasticMaterialStatus * ) this->giveStatus(gp);

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
    this->computeStressSpaceHardeningVars(fullStressSpaceHardeningVars, gp, strainSpaceHardeningVariables);

    //
    // check for plastic process
    //
    elastic = 1;
    actSurf = 0;
    activeConditionMap.zero();
    for ( i = 1; i <= this->nsurf; i++ ) {
        yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, fullStressSpaceHardeningVars);
        if ( yieldValue > YIELD_TOL ) {
            elastic = 0;
            actSurf += 1;
            activeConditionMap.at(i) = actSurf;
        }
    }

    if ( !elastic ) {
        /*
         * plastic multisurface closest-point corrector
         */


        do {
            do { // restart loop
                elasticStrainVectorR = totalStrain;
                elasticStrainVectorR.subtract(plasticStrainVectorR);
                // stress vector in full form due to computational convinience
                //if (fullStressVector) delete fullStressVector;
                this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
                this->computeStressSpaceHardeningVars(fullStressSpaceHardeningVars, gp, strainSpaceHardeningVariables);

                // compute gradients
                for ( i = 1; i <= nsurf; i++ ) {
                    this->computeGradientVector(yieldGradVec [ i - 1 ], yieldFunction, i, gp, fullStressVector, fullStressSpaceHardeningVars);
                }

                if ( this->plType == nonassociatedPT ) {
                    for ( i = 1; i <= nsurf; i++ ) {
                        this->computeGradientVector(loadGradVec [ i - 1 ], loadFunction, i, gp, fullStressVector, fullStressSpaceHardeningVars);
                    }
                }


                //this->computeGradientVector (gradientVectorR, gp, &fullStressVector, fullStressSpaceHardeningVars);
                this->computeResidualVector(residualVectorR, gp, gamma, activeConditionMap, plasticStrainVectorR,
                                            strainSpaceHardeningVariables, * loadGradVecPtr);

                // check convergence
                yieldConsistency = true;
                for ( i = 1; i <= nsurf; i++ ) {
                    if ( activeConditionMap.at(i) ) {
                        yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, fullStressSpaceHardeningVars);
                        if ( ( yieldValue >= YIELD_TOL ) ) {
                            yieldConsistency = false;
                        }
                    }
                }

                if ( yieldConsistency && ( residualVectorR.computeNorm() < RES_TOL ) ) {
                    answer = fullStressVector;
                    printf(" (%d iterations)", nIterations);
                    return;
                }

                // compute  consistent tangent moduli
                this->computeHardeningReducedModuli(hardeningModuli, gp, strainSpaceHardeningVariables, atTime);
                if ( hardeningModuli.giveNumberOfRows() ) {
                    hardeningModuliInverse.beInverseOf(hardeningModuli);
                    // delete hardeningModuli;
                } else {
                    hardeningModuliInverse.resize(0, 0);
                }

                this->computeAlgorithmicModuli(consistentModuli, gp, elasticModuliInverse,
                                               hardeningModuliInverse, gamma, activeConditionMap,
                                               fullStressVector, fullStressSpaceHardeningVars);

                rhs.resize(actSurf);
                gmat.resize(actSurf, actSurf);
                for ( i = 1; i <= nsurf; i++ ) {
                    if ( activeConditionMap.at(i) ) {
                        yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, fullStressSpaceHardeningVars);
                        helpVector.beTProductOf(consistentModuli, ( * yieldGradVecPtr ) [ i - 1 ]);

                        for ( j = 1; j <= this->nsurf; j++ ) {
                            if ( activeConditionMap.at(j) ) {
                                gmat.at(i, j) = ( * loadGradVecPtr ) [ j - 1 ].dotProduct(helpVector);
                            }
                        }

                        rhs.at(i) = yieldValue - residualVectorR.dotProduct( helpVector);
                    }
                }

                // obtain increment to consistency parameter
                gmat.solveForRhs(rhs, helpVector);

                // assign gammas
                for ( i = 1; i <= this->nsurf; i++ ) {
                    if ( ( indx = activeConditionMap.at(i) ) ) {
                        dgamma.at(i) = helpVector.at(indx);
                    } else {
                        dgamma.at(i) = 0.0;
                    }
                }

                tempGamma = gamma;
                tempGamma.add(dgamma);

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

                    continue;
                } else {
                    break;
                }
            } while ( 1 ); // end restart loop

            // obtain incremental plastic strains and internal variables
            // residualVectorR vector and gradVec array are changed
            // but they are computed once again when new iteration starts
            for ( i = 1; i <= nsurf; i++ ) {
                if ( activeConditionMap.at(i) ) {
                    ( * loadGradVecPtr ) [ i - 1 ].times( dgamma.at(i) );
                    residualVectorR.add( ( * loadGradVecPtr ) [ i - 1 ] );
                }
            }

            helpVector.beProductOf(consistentModuli, residualVectorR);
            this->computeDiagModuli(helpMtrx, gp, elasticModuliInverse, hardeningModuliInverse);
            helpVector2.beProductOf(helpMtrx, helpVector);

            // Update state variables and consistency parameter
            for ( i = 1; i <= strSize; i++ ) {
                plasticStrainVectorR.at(i) += helpVector2.at(i);
            }

            for ( i = strSize + 1; i <= totSize; i++ ) {
                strainSpaceHardeningVariables.at(i - strSize) += helpVector2.at(i);
            }

            gamma.add(dgamma);

            // increment iteration count
            nIterations++;

            if ( nIterations > PLASTIC_MATERIAL_MAX_ITERATIONS ) {
                _warning4( "GiveRealStressVector: local equlibrium not reached in %d iterations\nElement %d, gp %d, continuing",
                          PLASTIC_MATERIAL_MAX_ITERATIONS, gp->giveElement()->giveNumber(), gp->giveNumber() );
                answer = fullStressVector;
                break;
            }
        } while ( 1 );
    } else {
        answer = fullStressVector;
    }
}


void
MPlasticMaterial :: cuttingPlaneReturn(FloatArray &answer,
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
    FloatArray fullStressVector, fullStressSpaceHardeningVars;
    FloatArray fSigmaGradientVectorR, gGradientVectorR, helpVector, trialStressIncrement, rhs;
    FloatArray di, dj;
    FloatMatrix elasticModuli, hardeningModuli, dmat, helpMtrx, gmatInv, gradMat;
    std :: vector< FloatArray >yieldGradVec(this->nsurf), loadGradVec(this->nsurf), * yieldGradVecPtr, * loadGradVecPtr;
    FloatArray dgamma(this->nsurf);
    double yieldValue;
    int nIterations = 0;
    int size, sizeR, i, j, elastic, restart, actSurf, indx, iindx, jindx;

    MPlasticMaterialStatus *status = ( MPlasticMaterialStatus * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );

    if ( this->plType == associatedPT ) {
        yieldGradVecPtr = loadGradVecPtr = & yieldGradVec;
    } else {
        yieldGradVecPtr = & yieldGradVec;
        loadGradVecPtr  = & loadGradVec;
    }

    // huhu:
    status->givePlasticStrainVector(plasticStrainVectorR);
    status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);

    dgamma.resize(nsurf);
    dgamma.zero();
    gamma.resize(nsurf);
    gamma.zero();
    // compute elastic moduli
    this->computeReducedElasticModuli(elasticModuli, gp, atTime);
    // compute  hardening moduli and its inverse
    this->computeHardeningReducedModuli(hardeningModuli, gp, strainSpaceHardeningVariables, atTime);

    // initialize activeConditionMap
    elasticStrainVectorR = totalStrain;
    elasticStrainVectorR.subtract(plasticStrainVectorR);
    // stress vector in full form due to computational convinience
    this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
    crossSection->giveReducedCharacteristicVector(trialStressIncrement, gp, fullStressVector);
    trialStressIncrement.subtract( status->giveStressVector() );
    this->computeStressSpaceHardeningVars(fullStressSpaceHardeningVars, gp, strainSpaceHardeningVariables);

    elastic = 1;
    actSurf = 0;
    activeConditionMap.zero();
    for ( i = 1; i <= this->nsurf; i++ ) {
        yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, fullStressSpaceHardeningVars);
        if ( yieldValue > YIELD_TOL ) {
            elastic = 0;
            actSurf += 1;
            activeConditionMap.at(i) = actSurf;
        }
    }

    if ( !elastic ) {
        //
        // compute consistent moduli
        //
        sizeR = elasticModuli.giveNumberOfRows();
        size = sizeR + hardeningModuli.giveNumberOfRows();

        dmat.resize(size, size);
        dmat.zero();

        for ( i = 1; i <= sizeR; i++ ) {
            for ( j = 1; j <= sizeR; j++ ) {
                dmat.at(i, j) = elasticModuli.at(i, j);
            }
        }

        for ( i = sizeR + 1; i <= size; i++ ) {
            for ( j = sizeR + 1; j <= size; j++ ) {
                dmat.at(i, j) = hardeningModuli.at(i - sizeR, j - sizeR);
            }
        }

        do {            // local equilibrium loop
            // update  hardening moduli and its inverse
            this->computeHardeningReducedModuli(hardeningModuli, gp, strainSpaceHardeningVariables, atTime);
            for ( i = sizeR + 1; i <= size; i++ ) {
                for ( j = sizeR + 1; j <= size; j++ ) {
                    dmat.at(i, j) = hardeningModuli.at(i - sizeR, j - sizeR);
                }
            }

            do {              // loop for restart
                gmatInv.resize(actSurf, actSurf);
                gmatInv.zero();
                rhs.resize(actSurf);
                rhs.zero();

                // compute gradients
                for ( i = 1; i <= nsurf; i++ ) {
                    this->computeGradientVector(yieldGradVec [ i - 1 ], yieldFunction, i, gp, fullStressVector, fullStressSpaceHardeningVars);
                }

                if ( this->plType == nonassociatedPT ) {
                    for ( i = 1; i <= nsurf; i++ ) {
                        this->computeGradientVector(loadGradVec [ i - 1 ], loadFunction, i, gp, fullStressVector, fullStressSpaceHardeningVars);
                    }
                }



                for ( i = 1; i <= nsurf; i++ ) {
                    if ( ( iindx = activeConditionMap.at(i) ) ) {
                        rhs.at(iindx) = this->computeYieldValueAt(gp, i, fullStressVector, fullStressSpaceHardeningVars);
                        helpVector.beTProductOf(dmat, ( * yieldGradVecPtr ) [ i - 1 ]);

                        for ( j = 1; j <= this->nsurf; j++ ) {
                            if ( ( jindx = activeConditionMap.at(j) ) ) {
                                gmatInv.at(iindx, jindx) = ( * loadGradVecPtr ) [ j - 1 ].dotProduct(helpVector);
                            }
                        }
                    }
                }

                /*
                 *                              // compute plastic multipliers for active conds.
                 *      //computee gmatInv
                 *      gradMat.resize(actSurf, size);
                 *      for (i=1; i<=nsurf; i++) {
                 *        if ((iindx = activeConditionMap.at(i))) {
                 *          rhs.at (iindx) = this-> computeYieldValueAt (gp, i, fullStressVector, fullStressSpaceHardeningVars);
                 *          for (j=1; j<=size; j++) {
                 *            gradMat.at(iindx, j) = (*yieldGradVecPtr)[i-1].at(j);
                 *          }
                 *        }
                 *      }
                 *      helpMtrx.resize(actSurf, size);
                 *      helpMtrx.beProductOf (gradMat, dmat);
                 *      gmatInv.beProductTOf (helpMtrx, gradMat);
                 */
                // solve for plastic multiplier increments
                gmatInv.solveForRhs(rhs, helpVector);
                // assign gammas
                for ( i = 1; i <= this->nsurf; i++ ) {
                    if ( ( indx = activeConditionMap.at(i) ) ) {
                        dgamma.at(i) = helpVector.at(indx);
                    } else {
                        dgamma.at(i) = 0.0;
                    }
                }



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

                    continue;
                } else {
                    break;
                }
            } while ( 1 );           // end restart loop

            // compute plastic strain vector
            helpVector.resize(sizeR);
            helpVector.zero();
            for ( i = 1; i <= this->nsurf; i++ ) {
                if ( activeConditionMap.at(i) ) {
                    ( * loadGradVecPtr ) [ i - 1 ].times( dgamma.at(i) );
                    for ( j = 1; j <= sizeR; j++ ) {
                        helpVector.at(j) += ( * loadGradVecPtr ) [ i - 1 ].at(j);
                    }
                }
            }

            plasticStrainVectorR.add(helpVector);
            // update strain space hardening vector
            helpVector.resize(size - sizeR);
            helpVector.zero();
            for ( i = 1; i <= this->nsurf; i++ ) {
                if ( activeConditionMap.at(i) ) {
                    //gradVec[i-1].times(gamma.at(i));
                    for ( j = sizeR + 1; j <= size; j++ ) {
                        helpVector.at(j - sizeR) += ( * loadGradVecPtr ) [ i - 1 ].at(j);
                    }
                }
            }

            strainSpaceHardeningVariables.add(helpVector);
            // update plastic multipliers
            gamma.add(dgamma);

            elasticStrainVectorR = totalStrain;
            elasticStrainVectorR.subtract(plasticStrainVectorR);
            // stress vector in full form due to computational convinience
            this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
            this->computeStressSpaceHardeningVars(fullStressSpaceHardeningVars, gp, strainSpaceHardeningVariables);

            elastic = 1;
            for ( i = 1; i <= this->nsurf; i++ ) {
                yieldValue = this->computeYieldValueAt(gp, i, fullStressVector, fullStressSpaceHardeningVars);
                // check for end of iteration
                if ( yieldValue > YIELD_TOL ) {
                    elastic = 0;
                }
            }

            if ( elastic ) {
                break;
            }

            // increment iteration count
            nIterations++;

            if ( nIterations > PLASTIC_MATERIAL_MAX_ITERATIONS ) {
                char buff [ 200 ], buff1 [ 150 ];
                sprintf(buff, "GiveRealStressVector: local equlibrium not reached in %d iterations", PLASTIC_MATERIAL_MAX_ITERATIONS);
                sprintf( buff1, "Element %d, gp %d, continuing", gp->giveElement()->giveNumber(), gp->giveNumber() );
                _warning2(buff, buff1);
                //nIterations = 0; goto huhu;
                break;
            }
        } while ( 1 );
    }

    printf(" (%d iterations)", nIterations);
    answer = fullStressVector;
}


void
MPlasticMaterial :: computeGradientVector(FloatArray &answer, functType ftype, int isurf,
                                          GaussPoint *gp,
                                          const FloatArray &fullStressVector,
                                          const FloatArray &fullStressSpaceHardeningVars)
{
    /*
     * Computes gradient vector R in reduced form.
     * R = {\der{\sigma}{f_{n+1}}, \der{q}{f_{n+1}}}^T
     * Gradient vector contains partial derivatives of yield function
     * with respect to stresses and with respect to
     * strain space hardening variables
     *
     * Note: variablex with R posfix are in reduced stress-space.
     */
    FloatArray stressGradient, stressGradientR;
    FloatArray stressSpaceHardVarGradient;
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );
    int i, isize, size;

    this->computeStressGradientVector(stressGradient, ftype, isurf, gp, fullStressVector,
                                      fullStressSpaceHardeningVars);

    crossSection->giveReducedCharacteristicVector(stressGradientR, gp, stressGradient);

    this->computeStressSpaceHardeningVarsReducedGradient(stressSpaceHardVarGradient, ftype, isurf, gp,
                                                         fullStressVector, fullStressSpaceHardeningVars);

    isize = stressGradientR.giveSize();
    if ( stressSpaceHardVarGradient.isNotEmpty() ) {
        size = isize + stressSpaceHardVarGradient.giveSize();
    } else {
        size = isize;
    }

    answer.resize(size);
    for ( i = 1; i <= isize; i++ ) {
        answer.at(i) = stressGradientR.at(i);
    }

    for ( i = isize + 1; i <= size; i++ ) {
        answer.at(i) = stressSpaceHardVarGradient.at(i - isize);
    }
}


void
MPlasticMaterial :: computeResidualVector(FloatArray &answer, GaussPoint *gp, const FloatArray &gamma,
                                          const IntArray &activeConditionMap, const FloatArray &plasticStrainVectorR,
                                          const FloatArray &strainSpaceHardeningVariables,
                                          std :: vector< FloatArray > &gradientVectorR)
{
    /* Computes Residual vector for closes point projection algorithm */

    FloatArray oldPlasticStrainVectorR, oldStrainSpaceHardeningVariables;
    int i, j, isize, size;
    MPlasticMaterialStatus *status = ( MPlasticMaterialStatus * ) this->giveStatus(gp);

    isize = plasticStrainVectorR.giveSize();
    size =  gradientVectorR [ 0 ].giveSize();

    answer.resize(size);
    status->givePlasticStrainVector(oldPlasticStrainVectorR);
    status->giveStrainSpaceHardeningVars(oldStrainSpaceHardeningVariables);

    for ( i = 1; i <= isize; i++ ) {
        answer.at(i) = oldPlasticStrainVectorR.at(i) - plasticStrainVectorR.at(i);
    }

    for ( i = isize + 1; i <= size; i++ ) {
        answer.at(i) = oldStrainSpaceHardeningVariables.at(i - isize) - strainSpaceHardeningVariables.at(i - isize);
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
MPlasticMaterial :: computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
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
MPlasticMaterial :: computeAlgorithmicModuli(FloatMatrix &answer,
                                             GaussPoint *gp,
                                             const FloatMatrix &elasticModuliInverse,
                                             const FloatMatrix &hardeningModuliInverse,
                                             const FloatArray &gamma,
                                             const IntArray &activeConditionMap,
                                             const FloatArray &fullStressVector,
                                             const FloatArray &fullStressSpaceHardeningVars)
{
    /* returns consistent moduli in reduced form.
     * Note: elasticModuli and hardeningModuli will be inverted
     *
     * stressVector in full form
     */
    FloatMatrix gradientMatrix;
    FloatMatrix helpInverse;
    int i, j, isize, size;

    isize = elasticModuliInverse.giveNumberOfRows();
    if ( this->hasHardening() ) {
        size = isize + hardeningModuliInverse.giveNumberOfRows();
    } else {
        size = isize;
    }

    // assemble consistent moduli
    helpInverse.resize(size, size);
    helpInverse.zero();
    for ( i = 1; i <= nsurf; i++ ) {
        if ( activeConditionMap.at(i) ) {
            // ask gradientMatrix
            this->computeReducedGradientMatrix(gradientMatrix,  i, gp, fullStressVector,
                                               fullStressSpaceHardeningVars);

            gradientMatrix.times( gamma.at(i) );
            helpInverse.add(gradientMatrix);
        }
    }

    for ( i = 1; i <= isize; i++ ) {
        for ( j = 1; j <= isize; j++ ) {
            helpInverse.at(i, j) += elasticModuliInverse.at(i, j);
        }
    }

    if ( hardeningModuliInverse.giveNumberOfRows() > 0 ) {
        for ( i = isize + 1; i <= size; i++ ) {
            for ( j = isize + 1; j <= size; j++ ) {
                helpInverse.at(i, j) += hardeningModuliInverse.at(i - isize, j - isize);
            }
        }
    }

    if ( hasHardening() && ( hardeningModuliInverse.giveNumberOfRows() == 0 ) ) {
        FloatMatrix help;
        help.beInverseOf(helpInverse);
        answer.resize( isize + fullStressSpaceHardeningVars.giveSize(), isize + fullStressSpaceHardeningVars.giveSize() );
        answer.zero();
        for ( i = 1; i <= size; i++ ) {
            for ( j = 1; j <= size; j++ ) {
                answer.at(i, j) = help.at(i, j);
            }
        }
    } else {
        answer.beInverseOf(helpInverse);
    }
}

// ----------------------------------------------------------------------------//


void
MPlasticMaterial :: giveConsistentStiffnessMatrix(FloatMatrix &answer,
                                                  MatResponseForm form,
                                                  MatResponseMode mode,
                                                  GaussPoint *gp,
                                                  TimeStep *atTime)
{
    //
    // returns receiver material matrix for given reached state
    // described by temp variables.
    //

    FloatMatrix consistentModuli, elasticModuli, hardeningModuli;
    FloatMatrix consistentModuliInverse, elasticModuliInverse, hardeningModuliInverse;
    FloatMatrix consistentSubModuli, gradientMatrix, gmatInv, gmat, gradMat, helpMtrx, helpMtrx2, nmat, sbm1, answerR;
    FloatArray gradientVector, stressVector, fullStressVector;
    FloatArray stressSpaceHardeningVars;
    FloatArray strainSpaceHardeningVariables, helpVector;
    std :: vector< FloatArray >yieldGradVec(this->nsurf), loadGradVec(this->nsurf), * yieldGradVecPtr, * loadGradVecPtr;
    FloatArray helpVec;

    IntArray activeConditionMap, mask;
    FloatArray gamma;
    int size, sizeR, i, j, iindx, actSurf = 0;

    MPlasticMaterialStatus *status = ( MPlasticMaterialStatus * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );

    if ( this->plType == associatedPT ) {
        yieldGradVecPtr = loadGradVecPtr = & yieldGradVec;
    } else {
        yieldGradVecPtr = & yieldGradVec;
        loadGradVecPtr  = & loadGradVec;
    }

    // ask for plastic consistency parameter
    status->giveTempGamma(gamma);
    status->giveTempActiveConditionMap(activeConditionMap);
    //
    // check for elastic cases
    //
    if ( ( status->giveTempStateFlag() == MPlasticMaterialStatus :: PM_Elastic ) || ( status->giveTempStateFlag() == MPlasticMaterialStatus :: PM_Unloading ) ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, TangentStiffness,
                                                                    gp, atTime);
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
    sizeR = elasticModuliInverse.giveNumberOfRows();

    // compute  hardening moduli and its inverse
    this->computeHardeningReducedModuli(hardeningModuli, gp, strainSpaceHardeningVariables, atTime);
    if ( hardeningModuli.giveNumberOfRows() ) {
        hardeningModuliInverse.beInverseOf(hardeningModuli);
    } else {
        hardeningModuliInverse.resize(0, 0);
    }

    stressVector = status->giveStressVector();
    crossSection->giveFullCharacteristicVector(fullStressVector, gp, stressVector);
    status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
    this->computeStressSpaceHardeningVars(stressSpaceHardeningVars, gp, strainSpaceHardeningVariables);


    //
    // compute consistent moduli
    //
    this->computeAlgorithmicModuli(consistentModuli, gp, elasticModuliInverse, hardeningModuliInverse, gamma,
                                   activeConditionMap, fullStressVector, stressSpaceHardeningVars);

    //computee gmatInv
    for ( i = 1; i <= nsurf; i++ ) {
        this->computeGradientVector(yieldGradVec [ i - 1 ], yieldFunction, i, gp, fullStressVector, stressSpaceHardeningVars);
    }

    if ( this->plType == nonassociatedPT ) {
        for ( i = 1; i <= nsurf; i++ ) {
            this->computeGradientVector(loadGradVec [ i - 1 ], loadFunction, i, gp, fullStressVector, stressSpaceHardeningVars);
        }
    }


    size = consistentModuli.giveNumberOfColumns();
    gradMat.resize(actSurf, size);
    for ( i = 1; i <= nsurf; i++ ) {
        if ( ( iindx = activeConditionMap.at(i) ) ) {
            for ( j = 1; j <= size; j++ ) {
                gradMat.at(iindx, j) = ( * yieldGradVecPtr ) [ i - 1 ].at(j);
            }
        }
    }

    if ( this->plType == associatedPT ) {
        helpMtrx.resize(actSurf, size);
        helpMtrx.beProductTOf(consistentModuli, gradMat);
        gmatInv.beProductOf(gradMat, helpMtrx);
        gmat.beInverseOf(gmatInv);

        nmat.resize(sizeR, size);
        sbm1.beSubMatrixOf(consistentModuli, 1, sizeR, 1, size);
        nmat.beProductTOf(sbm1, gradMat);
        helpMtrx.beProductOf(nmat, gmat);
        sbm1.beSubMatrixOf(consistentModuli, 1, size, 1, sizeR);
        nmat.beProductOf(gradMat, sbm1);

        helpMtrx2.beProductOf(helpMtrx, nmat);

        answer.beSubMatrixOf(consistentModuli, 1, sizeR, 1, sizeR);
        answer.subtract(helpMtrx2);
    } else {
        FloatMatrix lgradMat(actSurf, size);
        for ( i = 1; i <= nsurf; i++ ) {
            if ( ( iindx = activeConditionMap.at(i) ) ) {
                for ( j = 1; j <= size; j++ ) {
                    lgradMat.at(iindx, j) = ( * loadGradVecPtr ) [ i - 1 ].at(j);
                }
            }
        }

        helpMtrx.resize(actSurf, size);
        helpMtrx.beProductTOf(consistentModuli, lgradMat);
        gmatInv.beProductOf(gradMat, helpMtrx);
        gmat.beInverseOf(gmatInv);

        nmat.resize(sizeR, size);
        sbm1.beSubMatrixOf(consistentModuli, 1, sizeR, 1, size);
        nmat.beProductTOf(sbm1, lgradMat);
        helpMtrx.beProductOf(nmat, gmat);
        sbm1.beSubMatrixOf(consistentModuli, 1, size, 1, sizeR);
        nmat.beProductOf(gradMat, sbm1);
        helpMtrx2.beProductOf(helpMtrx, nmat);

        answer.beSubMatrixOf(consistentModuli, 1, sizeR, 1, sizeR);
        answer.subtract(helpMtrx2);
    }

    if ( form == ReducedForm ) {
        return;
    } else {
        this->giveStressStrainMask( mask, ReducedForm, gp->giveMaterialMode() );
        answerR = answer;
        answer.beSubMatrixOfSizeOf(answerR, mask, 6);
    }
}


void
MPlasticMaterial :: giveElastoPlasticStiffnessMatrix(FloatMatrix &answer,
                                                     MatResponseForm form,
                                                     MatResponseMode mode,
                                                     GaussPoint *gp,
                                                     TimeStep *atTime)
{
    //
    // returns receiver material matrix for given reached state
    // described by temp variables.
    //

    FloatMatrix dmat, elasticModuli, hardeningModuli;
    FloatMatrix gradientMatrix, gmatInv, gmat, gradMat, helpMtrx, helpMtrx2;
    FloatArray gradientVector, stressVector, fullStressVector;
    FloatArray stressSpaceHardeningVars;
    FloatArray strainSpaceHardeningVariables, helpVector;
    std :: vector< FloatArray >yieldGradVec(this->nsurf), loadGradVec(this->nsurf), * yieldGradVecPtr, * loadGradVecPtr;
    FloatArray helpVec;

    IntArray activeConditionMap, mask;
    FloatArray gamma;
    int size, sizeR, i, j, iindx, actSurf = 0;

    MPlasticMaterialStatus *status = ( MPlasticMaterialStatus * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );

    if ( this->plType == associatedPT ) {
        yieldGradVecPtr = loadGradVecPtr = & yieldGradVec;
    } else {
        yieldGradVecPtr = & yieldGradVec;
        loadGradVecPtr  = & loadGradVec;
    }

    // ask for plastic consistency parameter
    status->giveTempGamma(gamma);
    status->giveTempActiveConditionMap(activeConditionMap);
    //
    // check for elastic cases
    //
    if ( ( status->giveTempStateFlag() == MPlasticMaterialStatus :: PM_Elastic ) || ( status->giveTempStateFlag() == MPlasticMaterialStatus :: PM_Unloading ) ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, TangentStiffness,
                                                                    gp, atTime);
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
    sizeR = elasticModuli.giveNumberOfRows();

    // compute  hardening moduli and its inverse
    this->computeHardeningReducedModuli(hardeningModuli, gp, strainSpaceHardeningVariables, atTime);

    stressVector = status->giveStressVector();
    crossSection->giveFullCharacteristicVector(fullStressVector, gp, stressVector);
    status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
    this->computeStressSpaceHardeningVars(stressSpaceHardeningVars, gp, strainSpaceHardeningVariables);


    //
    // compute consistent moduli
    //
    size = elasticModuli.giveNumberOfRows() + hardeningModuli.giveNumberOfRows();
    dmat.resize(size, size);
    dmat.zero();

    for ( i = 1; i <= sizeR; i++ ) {
        for ( j = 1; j <= sizeR; j++ ) {
            dmat.at(i, j) = elasticModuli.at(i, j);
        }
    }

    for ( i = sizeR + 1; i <= size; i++ ) {
        for ( j = sizeR + 1; j <= size; j++ ) {
            dmat.at(i, j) = hardeningModuli.at(i - sizeR, j - sizeR);
        }
    }


    //computee gmatInv
    for ( i = 1; i <= nsurf; i++ ) {
        this->computeGradientVector(yieldGradVec [ i - 1 ], yieldFunction, i, gp, fullStressVector, stressSpaceHardeningVars);
    }

    if ( this->plType == nonassociatedPT ) {
        for ( i = 1; i <= nsurf; i++ ) {
            this->computeGradientVector(loadGradVec [ i - 1 ], loadFunction, i, gp, fullStressVector, stressSpaceHardeningVars);
        }
    }

    gradMat.resize(actSurf, size);
    for ( i = 1; i <= nsurf; i++ ) {
        if ( ( iindx = activeConditionMap.at(i) ) ) {
            for ( j = 1; j <= size; j++ ) {
                gradMat.at(iindx, j) = yieldGradVec [ i - 1 ].at(j);
            }
        }
    }

    if ( this->plType == associatedPT ) {
        helpMtrx.resize(actSurf, size);
        helpMtrx.beProductOf(gradMat, dmat);
        gmatInv.beProductTOf(helpMtrx, gradMat);
        gmat.beInverseOf(gmatInv);

        helpMtrx.beProductOf(gradMat, elasticModuli);
        helpMtrx2.beProductOf(gmat, helpMtrx);
        answer.beTProductOf(helpMtrx, helpMtrx2);
        answer.negated();
        answer.add(elasticModuli);
    } else {
        FloatMatrix lgradMat(actSurf, size);
        for ( i = 1; i <= nsurf; i++ ) {
            if ( ( iindx = activeConditionMap.at(i) ) ) {
                for ( j = 1; j <= size; j++ ) {
                    lgradMat.at(iindx, j) = loadGradVec [ i - 1 ].at(j);
                }
            }
        }

        helpMtrx.resize(actSurf, size);
        helpMtrx.beProductOf(gradMat, dmat);
        gmatInv.beProductTOf(helpMtrx, lgradMat);
        gmat.beInverseOf(gmatInv);

        helpMtrx.beProductOf(gradMat, elasticModuli);
        helpMtrx2.beProductOf(gmat, helpMtrx);
        helpMtrx.beProductTOf(elasticModuli, lgradMat);
        answer.beTProductOf(helpMtrx, helpMtrx2);
        answer.negated();
        answer.add(elasticModuli);
    }

    if ( form == ReducedForm ) {
        return;
    } else {
        this->giveStressStrainMask( mask, ReducedForm, gp->giveMaterialMode() );
        helpMtrx = answer;
        answer.beSubMatrixOfSizeOf(helpMtrx, mask, 6);
    }
}


void
MPlasticMaterial :: computeDiagModuli(FloatMatrix &answer,
                                      GaussPoint *gp, FloatMatrix &elasticModuliInverse,
                                      FloatMatrix &hardeningModuliInverse)
{
    //
    // assembles diagonal moduli from elasticModuliInverse and hardeningModuliInverse
    //
    int size1, size2, i, j;

    size1 = elasticModuliInverse.giveNumberOfRows();
    if ( hardeningModuliInverse.giveNumberOfRows() ) {
        size2 = size1 + hardeningModuliInverse.giveNumberOfRows();
    } else {
        size2 = size1;
    }

    answer.resize(size2, size2);
    answer.zero();

    for ( i = 1; i <= size1; i++ ) {
        for ( j = 1; j <= size1; j++ ) {
            answer.at(i, j) = elasticModuliInverse.at(i, j);
        }
    }

    for ( i = size1 + 1; i <= size2; i++ ) {
        for ( j = size1 + 1; j <= size2; j++ ) {
            answer.at(i, j) = hardeningModuliInverse.at(i - size1, j - size1);
        }
    }
}


void
MPlasticMaterial :: computeReducedElasticModuli(FloatMatrix &answer,
                                                GaussPoint *gp,
                                                TimeStep *atTime)
{  /* Returns elastic moduli in reduced stress-strain space*/
    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, ReducedForm,
                                                                ElasticStiffness,
                                                                gp, atTime);
}


// overloaded from structural material

void
MPlasticMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form,
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
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    } else if ( rmType == mpm_ClosestPoint ) {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    } else {
        this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
MPlasticMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form,
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
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    } else if ( rmType == mpm_ClosestPoint ) {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    } else {
        this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
MPlasticMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                             MatResponseMode mode,
                                             GaussPoint *gp,
                                             TimeStep *atTime)

//
// return receiver's 2dPlaneStrainMtrx constructed from
// general 3dMatrialStiffnessMatrix
// (2dPlaneStrain ==> eps_z = gamma_xz = gamma_yz = 0.)
//
{
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    } else if ( rmType == mpm_ClosestPoint ) {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    } else {
        this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
MPlasticMaterial :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *atTime)

//
// returns receiver's 1dMaterialStiffnessMAtrix
// (1d case ==> sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.)
{
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    } else if ( rmType == mpm_ClosestPoint ) {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    } else {
        this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
MPlasticMaterial :: give2dBeamLayerStiffMtrx(FloatMatrix &answer, MatResponseForm form,
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
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    } else if ( rmType == mpm_ClosestPoint ) {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    } else {
        this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
MPlasticMaterial :: give2dPlateLayerStiffMtrx(FloatMatrix &answer,
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
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    } else if ( rmType == mpm_ClosestPoint ) {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    } else {
        this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
MPlasticMaterial :: give1dFiberStiffMtrx(FloatMatrix &answer,
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
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    } else if ( rmType == mpm_ClosestPoint ) {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    } else {
        this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
MPlasticMaterial :: give3dShellLayerStiffMtrx(FloatMatrix &answer, MatResponseForm form,
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
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    } else if ( rmType == mpm_ClosestPoint ) {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    } else {
        this->giveElastoPlasticStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


int
MPlasticMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    MPlasticMaterialStatus *status = ( MPlasticMaterialStatus * ) this->giveStatus(aGaussPoint);
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
MPlasticMaterial :: giveIPValueType(InternalStateType type)
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
MPlasticMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
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
MPlasticMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( type == IST_PlasticStrainTensor ) {
        return this->giveSizeOfReducedStressStrainVector( aGaussPoint->giveMaterialMode() );
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        return 3;
    } else {
        return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}


MPlasticMaterialStatus :: MPlasticMaterialStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), plasticStrainVector(), tempPlasticStrainVector(),
    strainSpaceHardeningVarsVector(), tempStrainSpaceHardeningVarsVector()
{
    state_flag = temp_state_flag = MPlasticMaterialStatus :: PM_Elastic;
    gamma = tempGamma = 0.;
}


MPlasticMaterialStatus :: ~MPlasticMaterialStatus()
{ }


void
MPlasticMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    int i, n;

    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( ( state_flag == MPlasticMaterialStatus :: PM_Yielding ) || ( state_flag == MPlasticMaterialStatus :: PM_Unloading ) ) {
        if ( state_flag == MPlasticMaterialStatus :: PM_Yielding ) {
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


void MPlasticMaterialStatus :: initTempStatus()
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
        strainSpaceHardeningVarsVector.resize( ( ( MPlasticMaterial * ) gp->giveMaterial() )->
                                              giveSizeOfReducedHardeningVarsVector(gp) );
        strainSpaceHardeningVarsVector.zero();
    }

    tempStrainSpaceHardeningVarsVector = strainSpaceHardeningVarsVector;

    temp_state_flag = state_flag;

    tempGamma = gamma;
    tempActiveConditionMap = activeConditionMap;
}


void
MPlasticMaterialStatus :: updateYourself(TimeStep *atTime)
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
MPlasticMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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
MPlasticMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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
