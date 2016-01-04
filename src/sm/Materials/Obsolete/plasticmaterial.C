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

#include "plasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "../sm/CrossSections/structuralcrosssection.h"
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
    return mode == _3dMat ||
           mode == _1dMat ||
           //<RESTRICTED_SECTION>
           mode == _PlaneStress ||
           //</RESTRICTED_SECTION>
           mode == _PlaneStrain ||
           mode == _PlateLayer ||
           mode == _2dBeamLayer ||
           mode == _Fiber;
}


MaterialStatus *
PlasticMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new PlasticMaterialStatus(1, this->giveDomain(), gp, this->giveSizeOfReducedHardeningVarsVector(gp));
}


void
PlasticMaterial :: giveRealStressVector(FloatArray &answer,
                                        GaussPoint *gp,
                                        const FloatArray &totalStrain,
                                        TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
// completely formulated in Reduced stress-strain space
{
    FloatArray strainSpaceHardeningVariables;
    FloatArray fullStressVector, *fullStressSpaceHardeningVars, *residualVectorR;
    FloatArray elasticStrainVectorR;
    FloatArray strainVectorR, plasticStrainVectorR, *gradientVectorR;
    FloatArray helpVec, helpVec2;
    double yieldValue, Gamma, dGamma, helpVal1, helpVal2;
    int strSize, totSize, nIterations = 0;
    FloatMatrix elasticModuli, hardeningModuli, consistentModuli;
    FloatMatrix elasticModuliInverse, hardeningModuliInverse;
    FloatMatrix helpMtrx, helpMtrx2;

    PlasticMaterialStatus *status = static_cast< PlasticMaterialStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(strainVectorR, gp, totalStrain,
                                                tStep, VM_Total);

    plasticStrainVectorR = status->givePlasticStrainVector();
    strainSpaceHardeningVariables = status->giveStrainSpaceHardeningVars();

    // tady konec debugovani - strainSpaceHardeningVariables ve statusu neinicializovany
    // to musi udelat material.

    Gamma = 0.;
    strSize = strainVectorR.giveSize(); // size of reducedStrain Vector
    totSize = strSize + strainSpaceHardeningVariables.giveSize();

    // compute elastic moduli and its inverse
    this->computeReducedElasticModuli(elasticModuli, gp, tStep);
    elasticModuliInverse.beInverseOf(elasticModuli);

    do {
        elasticStrainVectorR.beDifferenceOf(strainVectorR, plasticStrainVectorR);
        // stress vector in full form due to computational convenience
        this->computeTrialStressIncrement(fullStressVector, gp, elasticStrainVectorR, tStep);
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
        this->computeHardeningReducedModuli(hardeningModuli, gp, & strainSpaceHardeningVariables, tStep);
        if ( hardeningModuli.giveNumberOfRows() ) {
            hardeningModuliInverse.beInverseOf(hardeningModuli);
        } else {
            hardeningModuliInverse.clear();
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
        residualVectorR->add(* gradientVectorR);
        helpVec.beProductOf(consistentModuli, * residualVectorR);
        // note elasticModuli and hardeningModuli are yet inverted
        this->computeDiagModuli(helpMtrx, gp, elasticModuliInverse, hardeningModuliInverse);
        helpVec2.beProductOf(helpMtrx, helpVec);

        // Update state variables and consistency parameter
        for ( int i = 1; i <= strSize; i++ ) {
            plasticStrainVectorR.at(i) += helpVec2.at(i);
        }

        for ( int i = strSize + 1; i <= totSize; i++ ) {
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
            OOFEM_WARNING("local equlibrium not reached in %d iterations\nElement %d, gp %d, continuing",
                      PLASTIC_MATERIAL_MAX_ITERATIONS, gp->giveElement()->giveNumber(), gp->giveNumber() );
            break;
        }
    } while ( 1 );

    // update temp state variables in gp and associted material status

    status->letTempStrainVectorBe(totalStrain);
    StructuralMaterial :: giveReducedSymVectorForm( helpVec, fullStressVector, gp->giveMaterialMode() );
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

    answer = status->giveTempStressVector();
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
    int isize, size;

    stressGradient = this->ComputeStressGradient(gp, fullStressVector,
                                                 fullStressSpaceHardeningVars);

    StructuralMaterial :: giveReducedSymVectorForm( stressGradientR, * stressGradient, gp->giveMaterialMode() );

    stressSpaceHardVarGradient = this->ComputeStressSpaceHardeningVarsReducedGradient(gp, fullStressVector, fullStressSpaceHardeningVars);

    isize = stressGradientR.giveSize();
    if ( stressSpaceHardVarGradient ) {
        size = isize + stressSpaceHardVarGradient->giveSize();
    } else {
        size = isize;
    }

    answer = new FloatArray(size);
    for ( int i = 1; i <= isize; i++ ) {
        answer->at(i) = stressGradientR.at(i);
    }

    for ( int i = isize + 1; i <= size; i++ ) {
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
    int isize, size;
    PlasticMaterialStatus *status = static_cast< PlasticMaterialStatus * >( this->giveStatus(gp) );

    isize = plasticStrainVectorR->giveSize();
    size = gradientVectorR->giveSize();

    answer = new FloatArray(size);
    oldPlasticStrainVectorR = status->givePlasticStrainVector();
    oldStrainSpaceHardeningVariables = status->giveStrainSpaceHardeningVars();

    for ( int i = 1; i <= isize; i++ ) {
        answer->at(i) = oldPlasticStrainVectorR.at(i) - plasticStrainVectorR->at(i) + Gamma *gradientVectorR->at(i);
    }

    for ( int i = isize + 1; i <= size; i++ ) {
        answer->at(i) = oldStrainSpaceHardeningVariables.at(i - isize) - strainSpaceHardeningVariables->at(i - isize) + Gamma *gradientVectorR->at(i);
    }

    return answer;
}


void
PlasticMaterial :: computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                               const FloatArray &elasticStrainVectorR,
                                               TimeStep *tStep)
{
    /* Computes the full trial elastic stress vector */

    OOFEM_ERROR("Unable to compute trial stress increment");
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
    int isize, size;

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

    for ( int i = 1; i <= isize; i++ ) {
        for ( int j = 1; j <= isize; j++ ) {
            helpInverse.at(i, j) = elasticModuliInverse.at(i, j) + Gamma *gradientMatrix.at(i, j);
        }

        for ( int j = isize + 1; j <= size; j++ ) {
            helpInverse.at(i, j) = Gamma * gradientMatrix.at(i, j);
        }
    }

    for ( int i = isize + 1; i <= size; i++ ) {
        for ( int j = 1; j <= isize; j++ ) {
            helpInverse.at(i, j) = Gamma * gradientMatrix.at(i, j);
        }

        for ( int j = isize + 1; j <= size; j++ ) {
            helpInverse.at(i, j) = hardeningModuliInverse.at(i - isize, j - isize) + Gamma *gradientMatrix.at(i, j);
        }
    }

    answer.beInverseOf(helpInverse);
}

// ----------------------------------------------------------------------------//

void
PlasticMaterial :: giveConsistentStiffnessMatrix(FloatMatrix &answer,
                                                 MatResponseMode mode,
                                                 GaussPoint *gp,
                                                 TimeStep *tStep)
{
    //
    // returns receiver material matrix for given reached state
    // described by temp variables.
    //

    FloatMatrix consistentModuli, elasticModuli, hardeningModuli;
    FloatMatrix elasticModuliInverse, hardeningModuliInverse;
    FloatMatrix consistentSubModuli, answerR;
    FloatArray *gradientVector, stressVector, fullStressVector;
    FloatArray *stressSpaceHardeningVars;
    FloatArray strainSpaceHardeningVariables, helpVector;
    double s, Gamma;
    int sizeR;
    PlasticMaterialStatus *status = static_cast< PlasticMaterialStatus * >( this->giveStatus(gp) );

    // ask for plastic consistency parameter
    Gamma = status->giveTempPlasticConsistencyPrameter();

    // check for elastic cases
    if ( ( status->giveTempStateFlag() == PM_Elastic ) || ( status->giveTempStateFlag() == PM_Unloading ) ) {
        this->computeReducedElasticModuli(elasticModuli, gp, tStep);
        answer = elasticModuli;
        return;
    }

    //
    // plastic case
    //
    // compute elastic moduli and its inverse
    this->computeReducedElasticModuli(elasticModuli, gp, tStep);
    elasticModuliInverse.beInverseOf(elasticModuli);
    sizeR = elasticModuliInverse.giveNumberOfRows();

    // compute  hardening moduli and its inverse
    this->computeHardeningReducedModuli(hardeningModuli, gp, & strainSpaceHardeningVariables, tStep);
    if ( hardeningModuli.giveNumberOfRows() ) {
        hardeningModuliInverse.beInverseOf(hardeningModuli);
    } else {
        hardeningModuliInverse.clear();
    }

    stressVector = status->giveStressVector();
    StructuralMaterial :: giveFullSymVectorForm( fullStressVector, stressVector, gp->giveMaterialMode() );
    strainSpaceHardeningVariables = status->giveStrainSpaceHardeningVars();
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

    for ( int i = 1; i <= sizeR; i++ ) {
        for ( int j = 1; j <= sizeR; j++ ) {
            answerR.at(i, j) += ( 1. / s ) * helpVector.at(i) * helpVector.at(j);
        }
    }

    delete gradientVector;
    delete stressSpaceHardeningVars;

    answer =  answerR;
}

void
PlasticMaterial :: computeDiagModuli(FloatMatrix &answer,
                                     GaussPoint *gp, FloatMatrix &elasticModuliInverse,
                                     FloatMatrix &hardeningModuliInverse)
{
    //
    // assembles diagonal moduli from elasticModuliInverse and hardeningModuliInverse
    //
    int size1, size2;

    size1 = elasticModuliInverse.giveNumberOfRows();
    if ( hardeningModuliInverse.giveNumberOfRows() ) {
        size2 = size1 + hardeningModuliInverse.giveNumberOfRows();
    } else {
        size2 = size1;
    }

    answer.resize(size2, size2);
    answer.zero();

    for ( int i = 1; i <= size1; i++ ) {
        for ( int j = 1; j <= size1; j++ ) {
            answer.at(i, j) = elasticModuliInverse.at(i, j);
        }
    }

    for ( int i = size1 + 1; i <= size2; i++ ) {
        for ( int j = size1 + 1; j <= size2; j++ ) {
            answer.at(i, j) = hardeningModuliInverse.at(i - size1, j - size1);
        }
    }
}


void
PlasticMaterial :: computeReducedElasticModuli(FloatMatrix &answer,
                                               GaussPoint *gp,
                                               TimeStep *tStep)
{  /* Returns elastic moduli in reduced stress-strain space*/
    this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer,
                                                           ElasticStiffness,
                                                           gp, tStep);
}


// overloaded from structural material

void
PlasticMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                 MatResponseMode mode,
                                                 GaussPoint *gp,
                                                 TimeStep *tStep)
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
        OOFEM_ERROR("Different stressStrain mode encountered");
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
        this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
    } else {
        this->giveConsistentStiffnessMatrix(answer, mode, gp, tStep);
    }
}


void
PlasticMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer,
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
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    } else {
        this->giveConsistentStiffnessMatrix(answer, mode, gp, tStep);
    }
}


void
PlasticMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep)

//
// return receiver's 2dPlaneStrainMtrx constructed from
// general 3dMatrialStiffnessMatrix
// (2dPlaneStrain ==> eps_z = gamma_xz = gamma_yz = 0.)
//
{
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    } else {
        this->giveConsistentStiffnessMatrix(answer, mode, gp, tStep);
    }
}


void
PlasticMaterial :: give1dStressStiffMtrx(FloatMatrix &answer,
                                         MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *tStep)

//
// returns receiver's 1dMaterialStiffnessMAtrix
// (1d case ==> sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.)
{
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    } else {
        this->giveConsistentStiffnessMatrix(answer, mode, gp, tStep);
    }
}


void
PlasticMaterial :: give2dBeamLayerStiffMtrx(FloatMatrix &answer,
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
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    } else {
        this->giveConsistentStiffnessMatrix(answer, mode, gp, tStep);
    }
}


void
PlasticMaterial :: givePlateLayerStiffMtrx(FloatMatrix &answer,
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
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    } else {
        this->giveConsistentStiffnessMatrix(answer, mode, gp, tStep);
    }
}


void
PlasticMaterial :: giveFiberStiffMtrx(FloatMatrix &answer,
                                      MatResponseMode mode,
                                      GaussPoint *gp,
                                      TimeStep *tStep)
//
// returns receiver's Fiber
// (1dFiber ==> sigma_y = sigma_z = tau_yz = 0.)
//
// standard method from Material Class overloaded, because no inversion is needed.
// the reduction from 3d case will not work
// this implementation should be faster.
{
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    } else {
        this->giveConsistentStiffnessMatrix(answer, mode, gp, tStep);
    }
}


int
PlasticMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    PlasticMaterialStatus *status = static_cast< PlasticMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        const FloatArray ep = status->givePlasticStrainVector();
        ///@todo Fill in correct full form values here! This just adds zeros!
        StructuralMaterial :: giveFullSymVectorForm( answer, ep, gp->giveMaterialMode() );
        return 1;
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        FloatArray st(6);
        const FloatArray &s = status->givePlasticStrainVector();

        ///@todo Fill in correct full form values here! This just adds zeros!
        StructuralMaterial :: giveFullSymVectorForm( st, s, gp->giveMaterialMode() );

        this->computePrincipalValues(answer, st, principal_strain);
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}


PlasticMaterialStatus :: PlasticMaterialStatus(int n, Domain *d, GaussPoint *g, int statusSize) :
    StructuralMaterialStatus(n, d, g), plasticStrainVector(), tempPlasticStrainVector(),
    strainSpaceHardeningVarsVector(statusSize), tempStrainSpaceHardeningVarsVector(statusSize)
{
    state_flag = temp_state_flag = PM_Elastic;
    gamma = temp_gamma = 0.;
}


PlasticMaterialStatus :: ~PlasticMaterialStatus()
{ }


void
PlasticMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( ( state_flag == PM_Yielding ) || ( state_flag == PM_Unloading ) ) {
        if ( state_flag == PM_Yielding ) {
            fprintf(file, " Yielding, ");
        } else {
            fprintf(file, " Unloading, ");
        }

        fprintf(file, " plastic strains ");
        for ( auto &val : plasticStrainVector ) {
            fprintf( file, " %.4e", val );
        }

        if ( strainSpaceHardeningVarsVector.giveSize() ) {
            fprintf(file, ", strain space hardening vars ");
            for ( auto &val : strainSpaceHardeningVarsVector ) {
                fprintf( file, " %.4e", val );
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
        plasticStrainVector.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        plasticStrainVector.zero();
    }

    tempPlasticStrainVector = plasticStrainVector;

    tempStrainSpaceHardeningVarsVector = strainSpaceHardeningVarsVector;

    temp_state_flag = state_flag;
    temp_gamma = gamma;
}


void
PlasticMaterialStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    StructuralMaterialStatus :: updateYourself(tStep);

    plasticStrainVector = tempPlasticStrainVector;
    strainSpaceHardeningVarsVector = tempStrainSpaceHardeningVarsVector;

    state_flag = temp_state_flag;
    gamma = temp_gamma;
}


contextIOResultType
PlasticMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
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
    if ( ( iores = plasticStrainVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = strainSpaceHardeningVarsVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(state_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(gamma) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


contextIOResultType
PlasticMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = plasticStrainVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = strainSpaceHardeningVarsVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(state_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(gamma) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK; // return success
}

void PlasticMaterialStatus :: copyStateVariables(const MaterialStatus &iStatus)
{
    StructuralMaterialStatus :: copyStateVariables(iStatus);

    MaterialStatus &tmpStat = const_cast< MaterialStatus & >(iStatus);
    const PlasticMaterialStatus &plastStatus = dynamic_cast< PlasticMaterialStatus & >(tmpStat);

    plasticStrainVector = plastStatus.givePlasticStrainVector();
    tempPlasticStrainVector = plastStatus.giveTempPlasticStrainVector();

    strainSpaceHardeningVarsVector = plastStatus.giveStrainSpaceHardeningVars();
    tempStrainSpaceHardeningVarsVector = plastStatus.givetempStrainSpaceHardeningVarsVector();

    state_flag = plastStatus.giveStateFlag();
    temp_state_flag = plastStatus.giveTempStateFlag();

    gamma = plastStatus.givePlasticConsistencyPrameter();
    temp_gamma = plastStatus.giveTempPlasticConsistencyPrameter();
}

void PlasticMaterialStatus :: addStateVariables(const MaterialStatus &iStatus)
{
    printf("Entering PlasticMaterialStatus :: copyAddVariables().\n");
}

void PlasticMaterialStatus :: printYourself()
{
    printf("I am a PlasticMaterialStatus. plasticStrainVector: \n");
    plasticStrainVector.printYourself();
}
} // end namespace oofem
