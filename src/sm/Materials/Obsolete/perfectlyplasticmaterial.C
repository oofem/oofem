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

#include "perfectlyplasticmaterial.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
int
PerfectlyPlasticMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    return mode == _3dMat ||
           //<RESTRICTED_SECTION>
           mode == _PlaneStress ||
           //</RESTRICTED_SECTION>
           mode == _PlaneStrain ||
           mode == _1dMat ||
           mode == _PlateLayer ||
           mode == _2dBeamLayer;
}


void
PerfectlyPlasticMaterial :: giveRealStressVector(FloatArray &answer,
                                                 GaussPoint *gp,
                                                 const FloatArray &totalStrain,
                                                 TimeStep *tStep)
//
// returns  stress vector (in full or reduced form )  of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
//
// may be good idea for nonlinear materials overload this function in order to
// capture strain - stress history for proper modelling receiver behaviour
// and in order to capture possible failure of material
//
// Note: formulated in full stress strain space
{
    FloatArray elasticStressIncrement;
    FloatArray workStress, *yieldStressGrad, workStress2, stressVector3d;
    FloatArray mStrainIncrement3d, mStressElasticIncrement3d, PlasticStrainVector3d;
    FloatArray dSigmaIncrement3d, *stressCorrection, *loadingStressGrad;
    FloatArray helpArray;
    FloatArray plasticStrainIncrement3d, strainIncrement, reducedStrain, reducedStrainIncrement;
    FloatArray statusFullStressVector, statusFullPlasticVector;
    FloatArray plasticStrainVector;
    PerfectlyPlasticMaterialStatus *status = static_cast< PerfectlyPlasticMaterialStatus * >( this->giveStatus(gp) );

    FloatMatrix dp;
    StructuralCrossSection *crossSection = static_cast< StructuralCrossSection * >
                                           ( gp->giveElement()->giveCrossSection() );
    double f0, f1, f2, help, dLambda, r, r1, m;

    // init temp variables (of YC,LC,Material) at the beginning of step
    this->initTempStatus(gp);
    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain,
                                                tStep, VM_Total);

    reducedStrainIncrement.beDifferenceOf( reducedStrain, status->giveStrainVector() );

    StructuralMaterial :: giveFullSymVectorForm( strainIncrement, reducedStrainIncrement, gp->giveMaterialMode() );
    plasticStrainVector = status->givePlasticStrainVector();
    StructuralMaterial :: giveFullSymVectorForm( statusFullPlasticVector, plasticStrainVector, gp->giveMaterialMode() );
    StructuralMaterial :: giveFullSymVectorForm( statusFullStressVector, status->giveStressVector(), gp->giveMaterialMode() );

    //
    // Note : formulated in full stress strain space
    //
    this->computeTrialStressIncrement(elasticStressIncrement, gp, strainIncrement, tStep);
    //
    // calculate deltaSigmaPlastic
    //
    workStress = statusFullStressVector;
    workStress.add(elasticStressIncrement);

    f0 = this->computeYCValueAt(gp, & statusFullStressVector,
                                & statusFullPlasticVector);

    f1 = this->computeYCValueAt(gp, & workStress,
                                & statusFullPlasticVector);

    PlasticStrainVector3d = statusFullPlasticVector;


    if ( f1 >= 0. ) { // loading surface violated by the elastic trial set
        if ( f0 < 0. ) { // previously in elastic area
            // element not yielding - set the print status
            status->setTempYieldFlag(0);
            // compute scaling factor
            r1 = -f0 / ( f1 - f0 ); // linear interpolation in f (Zienkiewicz, 1969)
            yieldStressGrad = this->GiveYCStressGradient(gp, & statusFullStressVector,
                                                         & statusFullPlasticVector);
            crossSection->imposeStressConstrainsOnGradient(gp, yieldStressGrad);

            help = 0.;
            for ( int i = 1; i <= 6; i++ ) {
                help += yieldStressGrad->at(i) * elasticStressIncrement.at(i);
            }

            delete yieldStressGrad;

            workStress2 = elasticStressIncrement;
            workStress2.times(r1);
            workStress2.add(statusFullStressVector);

            f2 = this->computeYCValueAt(gp, & workStress2,
                                        & statusFullPlasticVector);

            if ( fabs(help) > 1.e-6 ) {
                r = r1 - f2 / help; // improved value of r (Nayak & Zienkiewicz, 1972)
            } else {
                r = r1;
            }
        } else { // f0 should be zero
            r  = 0.;
        }

        stressVector3d =  elasticStressIncrement;
        stressVector3d.times(r);
        stressVector3d.add(statusFullStressVector);

        m = STRAIN_STEPS;
        mStrainIncrement3d = strainIncrement;
        mStrainIncrement3d.times( ( 1.0 - r ) / m );
        mStressElasticIncrement3d = elasticStressIncrement;
        mStressElasticIncrement3d.times( ( 1.0 - r ) / m );

        // test if fracture or failure occurs
        this->updateIfFailure(gp, & stressVector3d, & PlasticStrainVector3d);
        // element yielding - set the print status
        status->setTempYieldFlag(1);
        // loop over m-steps
        for ( int i = 1; i <= m; i++ ) {
            //   yieldStressGrad = yieldCriteria->
            //    GiveStressGradient (gp, stressVector3d,
            //            PlasticStrainVector3d,
            //            gp-> giveHardeningParam());

            // compute dLambda and dp
            this->computePlasticStiffnessAt(dp, gp,
                                            & stressVector3d,
                                            & PlasticStrainVector3d,
                                            & mStrainIncrement3d,
                                            tStep,
                                            dLambda);
            if ( dLambda < 0. ) {
                dLambda = 0.;
            }

            dSigmaIncrement3d.beProductOf(dp, mStrainIncrement3d);
            dSigmaIncrement3d.negated();
            dSigmaIncrement3d.add(mStressElasticIncrement3d);
            // compute normal to loading surface
            loadingStressGrad = this->GiveLCStressGradient(gp, & stressVector3d,
                                                           & PlasticStrainVector3d);

            // update the stress state
            stressVector3d.add(dSigmaIncrement3d);
            // compute stress corrections

            stressCorrection = this->
                               GiveStressCorrectionBackToYieldSurface(gp, & stressVector3d,
                                                                      & PlasticStrainVector3d);
            // apply the stress correction -> correction is in the direction of
            // the normal to the yield surface
            stressVector3d.add(* stressCorrection);
            //   test = yieldCriteria-> computeValueAt(stressVector3d,
            //                 PlasticStrainVector3d,
            //                 gp-> giveHardeningParam());

            // calculate plastic strain increment

            crossSection->imposeStrainConstrainsOnGradient(gp, loadingStressGrad);
            loadingStressGrad->times(dLambda);
            PlasticStrainVector3d.add(* loadingStressGrad);
            // strainVector3d -> add (mStrainIncrement3d);

            // update yieldCriteria and loading criteria records to newly reached state
            // in loadingStressGrad is stored current plastic vector increment
            this->updateTempYC(gp, & stressVector3d, & PlasticStrainVector3d);
            this->updateTempLC(gp, & stressVector3d, & PlasticStrainVector3d);

            // test if fracture or failure occurs
            this->updateIfFailure(gp, & stressVector3d, & PlasticStrainVector3d);

            delete loadingStressGrad;
            delete stressCorrection;
        }

        // compute total stress increment during strainIncrement
        //  totalStressIncrement = statusFullStressVector;
        //  totalStressIncrement.times (-1.0);
        //  totalStressIncrement.add (stressVector3d);
        // update gp
        FloatArray helpArry;
        StructuralMaterial :: giveReducedSymVectorForm( helpArry, stressVector3d, gp->giveMaterialMode() );

        status->letTempStressVectorBe(helpArry);
        status->letTempStrainVectorBe(totalStrain);
    } else {
        // element not yielding - set the print status
        status->setTempYieldFlag(0);
        stressVector3d = statusFullStressVector;
        stressVector3d.add(elasticStressIncrement);
        // test if fracture or failure occurs
        this->updateIfFailure(gp, & stressVector3d, & PlasticStrainVector3d);
        // update gp
        status->letTempStrainVectorBe(totalStrain);
        StructuralMaterial :: giveReducedSymVectorForm( helpArray, stressVector3d, gp->giveMaterialMode() );
        status->letTempStressVectorBe(helpArray);
    }

    // update gp plastic strain
    plasticStrainIncrement3d.beDifferenceOf(PlasticStrainVector3d, statusFullPlasticVector);
    StructuralMaterial :: giveReducedSymVectorForm( helpArray, plasticStrainIncrement3d, gp->giveMaterialMode() );
    status->letPlasticStrainIncrementVectorBe(helpArray);

    StructuralMaterial :: giveReducedSymVectorForm( answer, stressVector3d, gp->giveMaterialMode() );
}


void
PerfectlyPlasticMaterial :: giveEffectiveMaterialStiffnessMatrix(FloatMatrix &answer,
                                                                 MatResponseMode mode,
                                                                 GaussPoint *gp,
                                                                 TimeStep *tStep)
//
//
// for case of perfectly plastic material
// computes full elastic constitutive matrix for case of gp stress-strain state.
// if strainIncrement == NULL a loading is assumed
//
// we follow terminology based on paper from R. de Borst:
// "Smeared Cracking, plasticity, creep - Unified Aproach"
//
// if derived material would like to implement failure behaviour
// it must redefine basic Give3dMaterialStiffnessMatrix function
// in order to take possible failure (tension cracking) into account
//
{
    // FloatMatrix *de; // elastic matrix respecting fracture or failure
    StructuralMaterial *lMat = static_cast< StructuralMaterial * >( this->giveLinearElasticMaterial() );

    if ( lMat->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
        FloatMatrix stiff;
        lMat->giveStiffnessMatrix(stiff, mode, gp, tStep);
        this->giveFullSymMatrixForm( answer, stiff, gp->giveMaterialMode() );
    } else {
        OOFEM_ERROR("giveEffectiveMaterialStiffnessMatrix - unsupported material mode");
    }
}

void
PerfectlyPlasticMaterial :: giveMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode,
                                                        GaussPoint *gp,
                                                        TimeStep *tStep)
//
//
//
// computes full constitutive matrix for case of gp stress-strain state.
// it returns elasto-plastic stiffness material matrix.
// if strainIncrement == NULL a loading is assumed
// for detailed description see (W.F.Chen: Plasticity in Reinforced Concrete, McGraw-Hill, 1982,
// chapter 6.)
//
// if derived material would like to implement failure behaviour
// it must redefine basic Give3dMaterialStiffnessMatrix function
// in order to take possible failure (tension cracking) into account
//
// in this function answer is Material stiffness matrix respecting
// current stress-strain mode in gp. This is reached by using
// impose constraints functions
{
    FloatArray statusFullStressVector, statusFullPlasticVector, plasticStrainVector;
    double lambda;
    PerfectlyPlasticMaterialStatus *status = static_cast< PerfectlyPlasticMaterialStatus * >( this->giveStatus(gp) );

    // double f = domain->giveYieldCriteria(yieldCriteria)->
    //  computeValueAt(gp->giveStressVector(), gp->givePlasticStrainVector(),
    //       gp-> giveHardeningParam());

    this->giveEffectiveMaterialStiffnessMatrix(answer, mode, gp, tStep);
    // if (f < YIELD_BOUNDARY)  // linear elastic behaviour
    //  return de;

    if ( status->giveYieldFlag() == 0 ) {
        return;
    }

    // if secant stiffness requested assume initial elastic matrix
    if ( mode == SecantStiffness ) {
        return;
    }

    plasticStrainVector = status->givePlasticStrainVector();
    StructuralMaterial :: giveFullSymVectorForm( statusFullPlasticVector, plasticStrainVector, gp->giveMaterialMode() );
    StructuralMaterial :: giveFullSymVectorForm( statusFullStressVector, status->giveStressVector(), gp->giveMaterialMode() );

    // yield condition satisfied
    FloatMatrix dp;
    this->computePlasticStiffnessAt(dp, gp, & statusFullStressVector,
                                    & statusFullPlasticVector,
                                    NULL,
                                    tStep,
                                    lambda);
    answer.add(dp);
}


void
PerfectlyPlasticMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
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
        this->giveMaterialStiffnessMatrix(answer, mode, gp, tStep);
    }
}


void
PerfectlyPlasticMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer,
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
    FloatMatrix fullAnswer;
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    } else {
        this->giveMaterialStiffnessMatrix(fullAnswer, mode, gp, tStep);
        StructuralMaterial :: giveReducedSymMatrixForm( answer, fullAnswer, gp->giveMaterialMode() );
    }
}


void
PerfectlyPlasticMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                                     MatResponseMode mode,
                                                     GaussPoint *gp,
                                                     TimeStep *tStep)

//
// return receiver's 2dPlaneStrainMtrx constructed from
// general 3dMatrialStiffnessMatrix
// (2dPlaneStrain ==> eps_z = gamma_xz = gamma_yz = 0.)
//
{
    FloatMatrix fullAnswer;
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    } else {
        this->giveMaterialStiffnessMatrix(fullAnswer, mode, gp, tStep);
        StructuralMaterial :: giveReducedSymMatrixForm( answer, fullAnswer, gp->giveMaterialMode() );
    }
}


void
PerfectlyPlasticMaterial :: give1dStressStiffMtrx(FloatMatrix &answer,
                                                  MatResponseMode mode,
                                                  GaussPoint *gp,
                                                  TimeStep *tStep)

//
// returns receiver's 1dMaterialStiffnessMAtrix
// (1d case ==> sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.)
{
    FloatMatrix fullAnswer;
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    } else {
        this->giveMaterialStiffnessMatrix(fullAnswer, mode, gp, tStep);
        StructuralMaterial :: giveReducedSymMatrixForm( answer, fullAnswer, gp->giveMaterialMode() );
    }
}



void
PerfectlyPlasticMaterial :: give2dBeamLayerStiffMtrx(FloatMatrix &answer,
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
    FloatMatrix fullAnswer;
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    } else {
        this->giveMaterialStiffnessMatrix(fullAnswer, mode, gp, tStep);
        StructuralMaterial :: giveReducedSymMatrixForm( answer, fullAnswer, gp->giveMaterialMode() );
    }
}


void
PerfectlyPlasticMaterial :: givePlateLayerStiffMtrx(FloatMatrix &answer,
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
    FloatMatrix fullAnswer;
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    } else {
        this->giveMaterialStiffnessMatrix(fullAnswer, mode, gp, tStep);
        StructuralMaterial :: giveReducedSymMatrixForm( answer, fullAnswer, gp->giveMaterialMode() );
    }
}


void
PerfectlyPlasticMaterial :: computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                                        const FloatArray &strainIncrement,
                                                        TimeStep *tStep)
//
// computest the elastic stress increment
// from stressIncrement in full stress strain space
//
{
    FloatMatrix materialMatrix;

    if ( strainIncrement.giveSize() == 0 ) {
        answer.clear();
        return;
    }

    this->giveEffectiveMaterialStiffnessMatrix(materialMatrix, TangentStiffness, gp, tStep);
    answer.beProductOf(materialMatrix, strainIncrement);
}


void
PerfectlyPlasticMaterial :: computePlasticStiffnessAt(FloatMatrix &answer,
                                                      GaussPoint *gp,
                                                      FloatArray *currentStressVector,
                                                      FloatArray *currentPlasticStrainVector,
                                                      FloatArray *strainIncrement3d,
                                                      TimeStep *tStep,
                                                      double &lambda)
//
// Computes full form of  Plastic stiffness Matrix at given state.
// gp is used only and only for setting proper MaterialMode ()
// returns proportionality factor lambda also if strainIncrement3d != NULL
{
    StructuralCrossSection *crossSection = static_cast< StructuralCrossSection * >
                                           ( gp->giveElement()->giveCrossSection() );
    FloatMatrix de, *yeldStressGradMat, *loadingStressGradMat;
    FloatMatrix fsde, gsfsde;
    FloatArray *yeldStressGrad, *loadingStressGrad;
    FloatArray help;
    double denominator, nominator;
    //
    // force de to be elastic even if gp in plastic state
    // to do this, a flag in this class exist -> ForceElasticResponce
    // if this flag is set to nonzero, function this::Give3dMaterialStiffnessMatrix
    // will be forced to return only elasticPart of 3dMaterialStiffnessMatrix
    // This function also set the ForceElasticFlagResponce to zero.
    //
    this->giveEffectiveMaterialStiffnessMatrix(de, TangentStiffness, gp, tStep);

    yeldStressGrad = this->GiveYCStressGradient(gp, currentStressVector,
                                                currentPlasticStrainVector);
    crossSection->imposeStressConstrainsOnGradient(gp, yeldStressGrad);
    yeldStressGradMat = new FloatMatrix(*yeldStressGrad, 1); // transpose

    loadingStressGrad = this->GiveLCStressGradient(gp, currentStressVector,
                                                   currentPlasticStrainVector);

    crossSection->imposeStrainConstrainsOnGradient(gp, loadingStressGrad);
    loadingStressGradMat = new FloatMatrix(*yeldStressGrad);

    help.beProductOf(de, * loadingStressGrad);
    delete loadingStressGrad;

    fsde.beProductOf(* yeldStressGradMat, de);
    delete yeldStressGradMat;
    gsfsde.beProductOf(* loadingStressGradMat, fsde);
    delete loadingStressGradMat;

    denominator = 0.;
    for ( int i = 1; i <= 6; i++ ) {
        denominator += yeldStressGrad->at(i) * help.at(i);
    }

    answer.beProductOf(de, gsfsde);

    answer.times( -( 1. / denominator ) );

    if ( strainIncrement3d != NULL ) { // compute proportional factor lambda
        nominator = 0.;
        help.beProductOf(de, * strainIncrement3d);
        for ( int i = 1; i <= 6; i++ ) {
            nominator += yeldStressGrad->at(i) * help.at(i);
        }

        lambda = nominator / denominator;
    }

    delete yeldStressGrad;
}


FloatArray *
PerfectlyPlasticMaterial :: GiveStressCorrectionBackToYieldSurface(GaussPoint *gp,
                                                                   FloatArray *stressVector3d,
                                                                   FloatArray *plasticVector3d)
//
// returns the stress correction -> correction is in the direction of
// the normal to the yield surface
// in full stress strain space
{
    FloatArray *yeldStressGrad, *stressCorrection;
    StructuralCrossSection *crossSection = static_cast< StructuralCrossSection * >
                                           ( gp->giveElement()->giveCrossSection() );
    double f3, help;

    yeldStressGrad = this->GiveYCStressGradient(gp,
                                                stressVector3d,
                                                plasticVector3d);
    crossSection->imposeStressConstrainsOnGradient(gp, yeldStressGrad);

    f3 = this->computeYCValueAt(gp, stressVector3d, plasticVector3d);

    help = 0.;
    for ( int j = 1; j <= 6; j++ ) {
        help += yeldStressGrad->at(j) * yeldStressGrad->at(j);
    }

    stressCorrection = new FloatArray(*yeldStressGrad);
    stressCorrection->times(-f3 / help);


    delete yeldStressGrad;
    return stressCorrection;
}


IRResultType
PerfectlyPlasticMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result = Material :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;
    return this->giveLinearElasticMaterial()->initializeFrom(ir);
}


double
PerfectlyPlasticMaterial :: give(int aProperty, GaussPoint *gp)
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
{
    double value = 0.0;

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


MaterialStatus *
PerfectlyPlasticMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates new  material status  corresponding to this class
 */
{
    return new PerfectlyPlasticMaterialStatus(1, this->giveDomain(), gp);
}


int
PerfectlyPlasticMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    PerfectlyPlasticMaterialStatus *status = static_cast< PerfectlyPlasticMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        const FloatArray &ep = status->givePlasticStrainVector();
        ///@todo Fill in correct full form values here! This just adds zeros!
        StructuralMaterial :: giveFullSymVectorForm( answer, ep, gp->giveMaterialMode() );
        return 1;
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        FloatArray st;

        const FloatArray &s = status->givePlasticStrainVector();
        ///@todo Fill in correct full form values here! This just adds zeros!
        StructuralMaterial :: giveFullSymVectorForm( st, s, gp->giveMaterialMode() );

        this->computePrincipalValues(answer, st, principal_strain);
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}


//##################################################################################################


PerfectlyPlasticMaterialStatus :: PerfectlyPlasticMaterialStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g),  plasticStrainVector(), plasticStrainIncrementVector()
{
    temp_yield_flag = yield_flag = 0; // elastic case at beginning
}


PerfectlyPlasticMaterialStatus :: ~PerfectlyPlasticMaterialStatus()
{ }


contextIOResultType
PerfectlyPlasticMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
//
{
    contextIOResultType iores;

    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream.write(yield_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = plasticStrainVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // return result back
    return CIO_OK;
}


contextIOResultType
PerfectlyPlasticMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restore state variables from stream
//
{
    contextIOResultType iores;

    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream.read(yield_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = plasticStrainVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // return result back
    return CIO_OK;
}


void
PerfectlyPlasticMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { yield_flag %d}\n", yield_flag);
}


void
PerfectlyPlasticMaterialStatus :: initTempStatus()
//
// initialize record at the begining of new load step
//
{
    StructuralMaterialStatus :: initTempStatus();

    if ( plasticStrainVector.giveSize() == 0 ) {
        plasticStrainVector.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        plasticStrainVector.zero();
    }

    if ( plasticStrainIncrementVector.giveSize() == 0 ) {
        plasticStrainIncrementVector.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        plasticStrainIncrementVector.zero();
    } else {
        plasticStrainIncrementVector.zero();
    }

    temp_yield_flag = yield_flag;
}


void
PerfectlyPlasticMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);

    plasticStrainVector.add(plasticStrainIncrementVector);
    plasticStrainIncrementVector.zero();

    yield_flag = temp_yield_flag;
}
} // end namespace oofem
