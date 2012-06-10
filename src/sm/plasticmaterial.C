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

#include "plasticmaterial.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "structuralcrosssection.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
#define YIELD_TOL 1.e-5
#define RES_TOL   1.e-5
#define PLASTIC_MATERIAL_MAX_ITERATIONS 40

PlasticMaterial :: PlasticMaterial(int n, Domain *d)  : StructuralMaterial(n, d)
//
// constructor
//
{
    linearElasticMaterial = NULL;
}


PlasticMaterial :: ~PlasticMaterial()
//
// destructor
//
{
    delete linearElasticMaterial;
}


int
PlasticMaterial :: hasMaterialModeCapability(MaterialMode mode)
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
PlasticMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates new  material status  corresponding to this class
 */
{
    PlasticMaterialStatus *status;

    status = new PlasticMaterialStatus(1, this->giveDomain(), gp);
    return status;
}


void
PlasticMaterial :: giveRealStressVector(FloatArray &answer,
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
    FloatArray fullStressVector, *fullStressSpaceHardeningVars, *residualVectorR;
    FloatArray strainIncrement, elasticStrainVectorR;
    FloatArray strainVectorR, plasticStrainVectorR, *gradientVectorR;
    FloatArray helpVec, helpVec2;
    double yieldValue, Gamma, dGamma, helpVal1, helpVal2;
    int i, strSize, totSize, nIterations = 0;
    FloatMatrix elasticModuli, hardeningModuli, consistentModuli;
    FloatMatrix elasticModuliInverse, hardeningModuliInverse;
    FloatMatrix helpMtrx, helpMtrx2;

    PlasticMaterialStatus *status = ( PlasticMaterialStatus * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(strainVectorR, gp, totalStrain,
                                                atTime, VM_Total);

    status->givePlasticStrainVector(plasticStrainVectorR);
    status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);

    // tady konec debugovani - strainSpaceHardeningVariables ve statusu neinicializovany
    // to musi udelat material.

    dGamma = Gamma = 0.;
    strSize = strainVectorR.giveSize(); // size of reducedStrain Vector
    totSize = strSize + strainSpaceHardeningVariables.giveSize();

    // compute elastic moduli and its inverse
    this->computeReducedElasticModuli(elasticModuli, gp, atTime);
    elasticModuliInverse.beInverseOf(elasticModuli);

    do {
        elasticStrainVectorR.beDifferenceOf(strainVectorR, plasticStrainVectorR);
        // stress vector in full form due to computational convenience
        this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, atTime);
        fullStressSpaceHardeningVars = this->ComputeStressSpaceHardeningVars(gp, & strainSpaceHardeningVariables);
        yieldValue = this->computeYieldValueAt(gp, & fullStressVector, fullStressSpaceHardeningVars);
        gradientVectorR = this->ComputeGradientVector(gp, & fullStressVector, fullStressSpaceHardeningVars);
        residualVectorR = this->ComputeResidualVector(gp, Gamma, & plasticStrainVectorR,
                                                      & strainSpaceHardeningVariables,
                                                      gradientVectorR);

        // check for end of iteration
        if ( ( yieldValue < YIELD_TOL ) && ( residualVectorR->computeNorm() < RES_TOL ) ) {
            delete fullStressSpaceHardeningVars;
            delete gradientVectorR;
            delete residualVectorR;
            break;
        }

        // compute  consistent tangent moduli
        this->computeHardeningReducedModuli(hardeningModuli, gp, & strainSpaceHardeningVariables, atTime);
        if ( hardeningModuli.giveNumberOfRows() ) {
            hardeningModuliInverse.beInverseOf(hardeningModuli);
        } else {
            hardeningModuliInverse.resize(0, 0);
        }

        this->computeConsistentModuli(consistentModuli, gp, elasticModuliInverse,
                                      hardeningModuliInverse, Gamma,
                                      fullStressVector, * fullStressSpaceHardeningVars);

        // obtain increment to consistency parameter
        helpMtrx.initFromVector(* gradientVectorR, 1);
        helpMtrx2.beProductOf(helpMtrx, consistentModuli);
        helpVec.beProductOf(helpMtrx2, * gradientVectorR);
        helpVal1 = helpVec.at(1);
        helpVec.beProductOf(helpMtrx2, * residualVectorR);
        helpVal2 = helpVec.at(1);
        dGamma = ( yieldValue - helpVal2 ) / helpVal1;

        // obtain incremental plastic strains and internal variables
        // we overwrite residualVectorR and gradientVectorR vectors
        // but they are computed once again when iteration continues
        gradientVectorR->times(dGamma);
        residualVectorR->add(*gradientVectorR);
        helpVec.beProductOf(consistentModuli, * residualVectorR);
        // note elasticModuli and hardeningModuli are yet inverted
        this->computeDiagModuli(helpMtrx, gp, elasticModuliInverse, hardeningModuliInverse);
        helpVec2.beProductOf(helpMtrx, helpVec);

        // Update state variables and consistency parameter
        for ( i = 1; i <= strSize; i++ ) {
            plasticStrainVectorR.at(i) += helpVec2.at(i);
        }

        for ( i = strSize + 1; i <= totSize; i++ ) {
            strainSpaceHardeningVariables.at(i - strSize) += helpVec2.at(i);
        }

        Gamma += dGamma;

        // increment iteration count
        nIterations++;

        // free allocated memory inside loop
        delete fullStressSpaceHardeningVars;
        delete gradientVectorR;
        delete residualVectorR;

        if ( nIterations > PLASTIC_MATERIAL_MAX_ITERATIONS ) {
            _warning4( "GiveRealStressVector: local equlibrium not reached in %d iterations\nElement %d, gp %d, continuing",
                      PLASTIC_MATERIAL_MAX_ITERATIONS, gp->giveElement()->giveNumber(), gp->giveNumber() );
            break;
        }
    } while ( 1 );

    // update temp state variables in gp and associted material status

    status->letTempStrainVectorBe(totalStrain);
    crossSection->giveReducedCharacteristicVector(helpVec, gp, fullStressVector);
    status->letTempStressVectorBe(helpVec);
    status->letTempPlasticStrainVectorBe(plasticStrainVectorR);
    status->letTempStrainSpaceHardeningVarsVectorBe(strainSpaceHardeningVariables);
    // update plastic consistency parameter
    status->letTempPlasticConsistencyPrameterBe(Gamma);

    // update state flag
    int newState, state = status->giveStateFlag();
    if ( Gamma > 0. ) {
        newState = PM_Yielding;            // test if plastic loading occur
    }
    // no plastic loading - check for unloading
    else if ( ( state == PM_Yielding ) || ( state == PM_Unloading ) ) {
        newState = PM_Unloading;
    } else {
        newState = PM_Elastic;
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


FloatArray *
PlasticMaterial :: ComputeGradientVector(GaussPoint *gp,
                                         FloatArray *fullStressVector,
                                         FloatArray *fullStressSpaceHardeningVars)
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
    FloatArray *stressGradient, stressGradientR;
    FloatArray *stressSpaceHardVarGradient, *answer;
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );
    int i, isize, size;

    stressGradient = this->ComputeStressGradient(gp, fullStressVector,
                                                 fullStressSpaceHardeningVars);

    crossSection->giveReducedCharacteristicVector(stressGradientR, gp, * stressGradient);

    stressSpaceHardVarGradient = this->ComputeStressSpaceHardeningVarsReducedGradient(gp, fullStressVector, fullStressSpaceHardeningVars);

    isize = stressGradientR.giveSize();
    if ( stressSpaceHardVarGradient ) {
        size = isize + stressSpaceHardVarGradient->giveSize();
    } else {
        size = isize;
    }

    answer = new FloatArray(size);
    for ( i = 1; i <= isize; i++ ) {
        answer->at(i) = stressGradientR.at(i);
    }

    for ( i = isize + 1; i <= size; i++ ) {
        answer->at(i) = stressSpaceHardVarGradient->at(i - isize);
    }

    delete stressSpaceHardVarGradient;
    delete stressGradient;

    return answer;
}


FloatArray *
PlasticMaterial :: ComputeResidualVector(GaussPoint *gp, double Gamma,
                                         FloatArray *plasticStrainVectorR,
                                         FloatArray *strainSpaceHardeningVariables,
                                         FloatArray *gradientVectorR)
{
    /* Computes Residual vector for closes point projection algorithm */

    FloatArray oldPlasticStrainVectorR, oldStrainSpaceHardeningVariables;
    FloatArray *answer;
    int i, isize, size;
    PlasticMaterialStatus *status = ( PlasticMaterialStatus * ) this->giveStatus(gp);

    isize = plasticStrainVectorR->giveSize();
    size = gradientVectorR->giveSize();

    answer = new FloatArray(size);
    status->givePlasticStrainVector(oldPlasticStrainVectorR);
    status->giveStrainSpaceHardeningVars(oldStrainSpaceHardeningVariables);

    for ( i = 1; i <= isize; i++ ) {
        answer->at(i) = oldPlasticStrainVectorR.at(i) - plasticStrainVectorR->at(i) + Gamma *gradientVectorR->at(i);
    }

    for ( i = isize + 1; i <= size; i++ ) {
        answer->at(i) = oldStrainSpaceHardeningVariables.at(i - isize) - strainSpaceHardeningVariables->at(i - isize) + Gamma *gradientVectorR->at(i);
    }

    return answer;
}


void
PlasticMaterial :: computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                               const FloatArray &elasticStrainVectorR,
                                               TimeStep *atTime)
{
    /* Computes the full trial elastic stress vector */

    _error("Unable to compute trial stress increment");
}


void
PlasticMaterial :: computeConsistentModuli(FloatMatrix &answer,
                                           GaussPoint *gp,
                                           FloatMatrix &elasticModuliInverse,
                                           FloatMatrix &hardeningModuliInverse,
                                           double Gamma,
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

    // ask gradientMatrix
    this->computeReducedGradientMatrix(gradientMatrix,  gp, fullStressVector,
                                       fullStressSpaceHardeningVars);

    // assemble consistent moduli
    helpInverse.resize(size, size);

    for ( i = 1; i <= isize; i++ ) {
        for ( j = 1; j <= isize; j++ ) {
            helpInverse.at(i, j) = elasticModuliInverse.at(i, j) + Gamma *gradientMatrix.at(i, j);
        }

        for ( j = isize + 1; j <= size; j++ ) {
            helpInverse.at(i, j) = Gamma * gradientMatrix.at(i, j);
        }
    }

    for ( i = isize + 1; i <= size; i++ ) {
        for ( j = 1; j <= isize; j++ ) {
            helpInverse.at(i, j) = Gamma * gradientMatrix.at(i, j);
        }

        for ( j = isize + 1; j <= size; j++ ) {
            helpInverse.at(i, j) = hardeningModuliInverse.at(i - isize, j - isize) + Gamma *gradientMatrix.at(i, j);
        }
    }

    answer.beInverseOf(helpInverse);
}

// ----------------------------------------------------------------------------//

void
PlasticMaterial :: giveConsistentStiffnessMatrix(FloatMatrix &answer,
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
    FloatMatrix consistentSubModuli, answerR;
    FloatArray *gradientVector, stressVector, fullStressVector;
    FloatArray *stressSpaceHardeningVars;
    FloatArray strainSpaceHardeningVariables, helpVector;
    IntArray mask;
    double s, Gamma;
    int sizeR, i, j;
    PlasticMaterialStatus *status = ( PlasticMaterialStatus * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );

    // ask for plastic consistency parameter
    Gamma = status->giveTempPlasticConsistencyPrameter();

    // check for elastic cases
    if ( ( status->giveTempStateFlag() == PM_Elastic ) || ( status->giveTempStateFlag() == PM_Unloading ) ) {
        this->computeReducedElasticModuli(elasticModuli, gp, atTime);
        if ( form == ReducedForm ) {
            answer = elasticModuli;
            return;
        } else {
            this->giveStressStrainMask( mask, ReducedForm, gp->giveMaterialMode() );
            answer.beSubMatrixOfSizeOf(elasticModuli, mask, 6);
            return;
        }
    }

    //
    // plastic case
    //
    // compute elastic moduli and its inverse
    this->computeReducedElasticModuli(elasticModuli, gp, atTime);
    elasticModuliInverse.beInverseOf(elasticModuli);
    sizeR = elasticModuliInverse.giveNumberOfRows();

    // compute  hardening moduli and its inverse
    this->computeHardeningReducedModuli(hardeningModuli, gp, & strainSpaceHardeningVariables, atTime);
    if ( hardeningModuli.giveNumberOfRows() ) {
        hardeningModuliInverse.beInverseOf(hardeningModuli);
    } else {
        hardeningModuliInverse.resize(0, 0);
    }

    stressVector = status->giveStressVector();
    crossSection->giveFullCharacteristicVector(fullStressVector, gp, stressVector);
    status->giveStrainSpaceHardeningVars(strainSpaceHardeningVariables);
    stressSpaceHardeningVars = this->ComputeStressSpaceHardeningVars(gp, & strainSpaceHardeningVariables);

    this->computeConsistentModuli(consistentModuli, gp,
                                  elasticModuliInverse,
                                  hardeningModuliInverse,
                                  Gamma,
                                  fullStressVector,
                                  * stressSpaceHardeningVars);

    gradientVector = this->ComputeGradientVector(gp, & fullStressVector, stressSpaceHardeningVars);
    helpVector.beProductOf(consistentModuli, * gradientVector);
    s = ( -1. ) * gradientVector->dotProduct(helpVector);

    answerR.beSubMatrixOf(consistentModuli, 1, sizeR, 1, sizeR);
    //consistentSubModuli.beSubMatrixOf (consistentModuli, 1,gradientVector->giveSize(), 1, sizeR);
    consistentSubModuli.beSubMatrixOf( consistentModuli, 1, sizeR, 1, gradientVector->giveSize() );
    helpVector.beProductOf(consistentSubModuli, * gradientVector);

    for ( i = 1; i <= sizeR; i++ ) {
        for ( j = 1; j <= sizeR; j++ ) {
            answerR.at(i, j) += ( 1. / s ) * helpVector.at(i) * helpVector.at(j);
        }
    }

    delete gradientVector;
    delete stressSpaceHardeningVars;

    if ( form == ReducedForm ) {
        answer =  answerR;
    } else {
        this->giveStressStrainMask( mask, ReducedForm, gp->giveMaterialMode() );
        answer.beSubMatrixOfSizeOf(answerR, mask, 6);
    }
}

void
PlasticMaterial :: computeDiagModuli(FloatMatrix &answer,
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
PlasticMaterial :: computeReducedElasticModuli(FloatMatrix &answer,
                                               GaussPoint *gp,
                                               TimeStep *atTime)
{  /* Returns elastic moduli in reduced stress-strain space*/
    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, ReducedForm,
                                                                ElasticStiffness,
                                                                gp, atTime);
}


// overloaded from structural material

void
PlasticMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form,
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
    } else {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
PlasticMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form,
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
    } else {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
PlasticMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseForm form,
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
    } else {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
PlasticMaterial :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                         MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *atTime)

//
// returns receiver's 1dMaterialStiffnessMAtrix
// (1d case ==> sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.)
{
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    } else {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
PlasticMaterial :: give2dBeamLayerStiffMtrx(FloatMatrix &answer, MatResponseForm form,
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
    } else {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
PlasticMaterial :: give2dPlateLayerStiffMtrx(FloatMatrix &answer,
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
    } else {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
PlasticMaterial :: give1dFiberStiffMtrx(FloatMatrix &answer,
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
    } else {
        this->giveConsistentStiffnessMatrix(answer, form, mode, gp, atTime);
    }
}


void
PlasticMaterial :: give3dShellLayerStiffMtrx(FloatMatrix &answer, MatResponseForm form,
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
    } else {
        this->give2dPlateLayerStiffMtrx(answer, form, mode, gp, atTime);
    }
}


int
PlasticMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    PlasticMaterialStatus *status = ( PlasticMaterialStatus * ) this->giveStatus(aGaussPoint);
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
PlasticMaterial :: giveIPValueType(InternalStateType type)
{
    // strains components packed in engineering notation
    if ( ( type == IST_PlasticStrainTensor ) || ( type == IST_PrincipalPlasticStrainTensor ) ) {
        return ISVT_TENSOR_S3E;
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }
}


int
PlasticMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
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
PlasticMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( type == IST_PlasticStrainTensor ) {
        return this->giveSizeOfReducedStressStrainVector( aGaussPoint->giveMaterialMode() );
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        return 3;
    } else {
        return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}


PlasticMaterialStatus :: PlasticMaterialStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), plasticStrainVector(), tempPlasticStrainVector(),
    strainSpaceHardeningVarsVector(), tempStrainSpaceHardeningVarsVector()
{
    state_flag = temp_state_flag = PM_Elastic;
    gamma = temp_gamma = 0.;
}


PlasticMaterialStatus :: ~PlasticMaterialStatus()
{ }


void
PlasticMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    int i, n;

    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( ( state_flag == PM_Yielding ) || ( state_flag == PM_Unloading ) ) {
        if ( state_flag == PM_Yielding ) {
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


void PlasticMaterialStatus :: initTempStatus()
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
        strainSpaceHardeningVarsVector.resize( ( ( PlasticMaterial * ) gp->giveMaterial() )->
                                              giveSizeOfReducedHardeningVarsVector(gp) );
        strainSpaceHardeningVarsVector.zero();
    }

    tempStrainSpaceHardeningVarsVector = strainSpaceHardeningVarsVector;

    temp_state_flag = state_flag;
    temp_gamma = gamma;
}


void
PlasticMaterialStatus :: updateYourself(TimeStep *atTime)
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
    gamma = temp_gamma;
}


contextIOResultType
PlasticMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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

    if ( !stream->write(& gamma, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


contextIOResultType
PlasticMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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

    if ( !stream->read(& gamma, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK; // return success
}

} // end namespace oofem
