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
 *               Copyright (C) 1993 - 2019   Borek Patzak
 *
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

#include "latticeframesteelplastic.h"
#include "latticeframeelastic.h"
#include "latticestructuralmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "CrossSections/structuralcrosssection.h"
#include "engngm.h"
#include "mathfem.h"
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "datastream.h"
#include "staggeredproblem.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(LatticeFrameSteelPlastic);

// constructor which creates a dummy material without a status and without random extension interface
// LatticeFrameSteelPlastic :: LatticeFrameSteelPlastic(int n, Domain *d) :
// LatticeStructuralMaterial(n, d)
bool
LatticeFrameSteelPlastic::hasMaterialModeCapability(MaterialMode mode) const
{
    return ( mode == _3dLattice );
}


void
LatticeFrameSteelPlastic::initializeFrom(InputRecord &ir)
{
    LatticeStructuralMaterial::initializeFrom(ir);

    //Young's modulus of the material that the beam element is made of
    IR_GIVE_FIELD(ir, this->e, _IFT_LatticeFrameSteelPlastic_e); // Macro

    //Poisson's ratio of the material that the beam element is made of
    IR_GIVE_FIELD(ir, this->nu, _IFT_LatticeFrameSteelPlastic_n); // Macro

    //Nx0
    IR_GIVE_FIELD(ir, this->nx0, _IFT_LatticeFrameSteelPlastic_nx0); // Macro

    //Mx0
    IR_GIVE_FIELD(ir, this->mx0, _IFT_LatticeFrameSteelPlastic_mx0); // Macro

    //My0
    IR_GIVE_FIELD(ir, this->my0, _IFT_LatticeFrameSteelPlastic_my0); // Macro

    //Mz0
    IR_GIVE_FIELD(ir, this->mz0, _IFT_LatticeFrameSteelPlastic_mz0); // Macro

    yieldTol = 1.e-6;
    ;
    IR_GIVE_FIELD(ir, this->yieldTol, _IFT_LatticeFrameSteelPlastic_tol); // Macro

    this->newtonIter = 100;
    IR_GIVE_FIELD(ir, this->newtonIter, _IFT_LatticeFrameSteelPlastic_iter); // Macro

    numberOfSubIncrements = 10;
    IR_GIVE_FIELD(ir, this->numberOfSubIncrements, _IFT_LatticeFrameSteelPlastic_sub); // Macro

    this->plasticFlag = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, plasticFlag, _IFT_LatticeFrameSteelPlastic_plastic); // Macro
}

MaterialStatus *
LatticeFrameSteelPlastic::CreateStatus(GaussPoint *gp) const
{
    return new LatticeFrameSteelPlasticStatus(1, LatticeFrameSteelPlastic::domain, gp);
}

MaterialStatus *
LatticeFrameSteelPlastic::giveStatus(GaussPoint *gp) const
{
    // test
    MaterialStatus *status = static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    if ( !status ) {
        // create a new one
        status = this->CreateStatus(gp);

        if ( status ) {
            gp->setMaterialStatus(status);
        }
    }

    return status;
}

double
LatticeFrameSteelPlastic::computeYieldValue(const FloatArrayF< 4 > &stress,
                                            GaussPoint *gp,
                                            TimeStep *tStep) const
{
    double yieldValue = 0.;
    double nx = stress.at(1);
    double mx = stress.at(2);
    double my = stress.at(3);
    double mz = stress.at(4);

    {
        yieldValue = pow(nx / this->nx0, 2.) + pow(mx / this->mx0, 2.) + pow(my / this->my0, 2.) + pow(mz / this->mz0, 2.) - 1.;
    }

    return yieldValue;
}

FloatArrayF< 4 >
LatticeFrameSteelPlastic::computeFVector(const FloatArrayF< 4 > &stress,
                                         GaussPoint *gp,
                                         TimeStep *tStep) const
{
    double nx = stress.at(1);
    double mx = stress.at(2);
    double my = stress.at(3);
    double mz = stress.at(4);

    FloatArrayF< 4 >f;

    f.at(1) = 2. * nx / pow(this->nx0, 2.);
    f.at(2) = 2. * mx / pow(this->mx0, 2.);
    f.at(3) = 2. * my / pow(this->my0, 2.);
    f.at(4) = 2. * mz / pow(this->mz0, 2.);

    return f;
}

FloatMatrixF< 4, 4 >
LatticeFrameSteelPlastic::computeDMMatrix(const FloatArrayF< 4 > &stress, GaussPoint *gp, TimeStep *tStep) const
{
    FloatMatrixF< 4, 4 >dm;

    //Derivatives of dGDSig
    dm.at(1, 1) = 2. / pow(this->nx0, 2.);
    dm.at(1, 2) = 0;
    dm.at(1, 3) = 0;
    dm.at(1, 4) = 0;

    //Derivatives of dGDTau
    dm.at(2, 1) = 0;
    dm.at(2, 2) = 2. / pow(this->mx0, 2.);
    dm.at(2, 3) = 0;
    dm.at(2, 4) = 0;

    //Derivates of evolution law
    dm.at(3, 1) = 0;
    dm.at(3, 2) = 0;
    dm.at(3, 3) = 2. / pow(this->my0, 2.);
    dm.at(3, 4) = 0;

    //Derivates of evolution law
    dm.at(4, 1) = 0;
    dm.at(4, 2) = 0;
    dm.at(4, 3) = 0;
    dm.at(4, 3) = 2. / pow(this->mz0, 2.);

    return dm;
}

FloatArrayF< 6 >
LatticeFrameSteelPlastic::giveThermalDilatationVector(GaussPoint *gp,  TimeStep *tStep) const
// returns a FloatArray(6) of initial strain vector caused by unit temperature in direction of gp (element) local axes
{
    double alpha = this->give(tAlpha, gp);

    return {
        alpha, 0., 0., 0., 0., 0.
    };
}

FloatArrayF< 6 >
LatticeFrameSteelPlastic::giveReducedStrain(GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );
    return status->giveReducedLatticeStrain();
}


FloatArrayF< 6 >
LatticeFrameSteelPlastic::performPlasticityReturn(GaussPoint *gp, const FloatArrayF< 6 > &reducedStrain, TimeStep *tStep) const
{
    double g = this->e / ( 2. * ( 1. + this->nu ) );
    const double area = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveArea();
    const double iy = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIy();
    const double iz = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIz();
    const double ik = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIk();

    auto status = static_cast< LatticeFrameSteelPlasticStatus * >( this->giveStatus(gp) );

    //Shear components are not used for plasticity return
    auto strain = reducedStrain [ { 0, 3, 4, 5 } ];

    /* Get plastic strain vector from status*/
    auto tempPlasticStrain = status->givePlasticLatticeStrain() [ { 0, 3, 4, 5 } ];

    FloatArrayF< 4 >tangent = { area *this->e, ik *g, iy *this->e, iz *this->e };

    /* Compute trial stress*/
    auto stress = mult(tangent, strain - tempPlasticStrain);

    auto oldStrain = this->giveReducedStrain(gp, tStep) [ { 0, 3, 4, 5 } ];

    /* Compute yield value*/
    double yieldValue = computeYieldValue(stress, gp, tStep);
    int subIncrementCounter = 0;

    /* Check yield condition, i.e. if the yield value is less than the yield tolerance. If yield condition is valid. Do perform regular return (closest point return)*/

    if ( yieldValue > yieldTol ) {
        int subIncrementFlag = 0;
        auto convergedStrain = oldStrain;
        auto tempStrain = strain;
        auto deltaStrain = strain - oldStrain;
        //To get into the loop
        status->letTempReturnResultBe(LatticeFrameSteelPlastic::RR_NotConverged);
        while ( status->giveTempReturnResult() == RR_NotConverged || subIncrementFlag == 1 ) {
            stress = mult(tangent, tempStrain - tempPlasticStrain);
            performRegularReturn(stress, yieldValue, gp, tStep);

            if ( status->giveTempReturnResult() == RR_NotConverged ) {
                subIncrementCounter++;
                if ( subIncrementCounter > numberOfSubIncrements ) {
                    OOFEM_LOG_INFO("Unstable element %d \n", gp->giveElement()->giveGlobalNumber() );
                    OOFEM_LOG_INFO("Yield value %e \n", yieldValue);
                    OOFEM_LOG_INFO("ConvergedStrain value %e %e %e %e\n", convergedStrain.at(1), convergedStrain.at(2), convergedStrain.at(3), convergedStrain.at(4) );
                    OOFEM_LOG_INFO("tempStrain value %e %e %e %e\n", tempStrain.at(1), tempStrain.at(2), tempStrain.at(3), tempStrain.at(4) );
                    OOFEM_LOG_INFO("deltaStrain value %e %e %e %e\n", deltaStrain.at(1), deltaStrain.at(2), deltaStrain.at(3), deltaStrain.at(4) );
                    OOFEM_LOG_INFO("targetstrain value %e %e %e %e\n", strain.at(1), strain.at(2), strain.at(3), strain.at(4) );

                    OOFEM_ERROR("LatticeFrameSteelPlastic :: performPlasticityReturn - Could not reach convergence with small deltaStrain, giving up.");
                }
                subIncrementFlag = 1;
                deltaStrain *= 0.5;
                tempStrain = convergedStrain + deltaStrain;
            } else if ( status->giveTempReturnResult() == RR_Converged && subIncrementFlag == 1 ) {
                tempPlasticStrain.at(1) = tempStrain.at(1) - stress.at(1) / ( area * this->e );
                tempPlasticStrain.at(2) = tempStrain.at(2) - stress.at(2) / ( ik * g );
                tempPlasticStrain.at(3) = tempStrain.at(3) - stress.at(3) / ( iy * this->e );
                tempPlasticStrain.at(4) = tempStrain.at(4) - stress.at(4) / ( iz * this->e );

                status->letTempPlasticLatticeStrainBe(assemble< 6 >(tempPlasticStrain, { 0, 3, 4, 5 }) );

                subIncrementFlag = 0;

                status->letTempReturnResultBe(LatticeFrameSteelPlastic::RR_NotConverged);
                convergedStrain = tempStrain;
                deltaStrain = strain - convergedStrain;
                tempStrain = strain;
                subIncrementCounter = 0;
            }
        }
    }

    const double shearareay = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveShearAreaY();
    const double shearareaz = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveShearAreaZ();


    tempPlasticStrain.at(1) = strain.at(1) - stress.at(1) / ( area * this->e );
    tempPlasticStrain.at(2) = strain.at(2) - stress.at(2) / ( ik * g );
    tempPlasticStrain.at(3) = strain.at(3) - stress.at(3) / ( iy * this->e );
    tempPlasticStrain.at(4) = strain.at(4) - stress.at(4) / ( iz * this->e );


    status->letTempPlasticLatticeStrainBe(assemble< 6 >(tempPlasticStrain, { 0, 3, 4, 5 }) );
    auto answer = assemble< 6 >(stress, { 0, 3, 4, 5 });
    answer.at(2) = shearareay * g * reducedStrain.at(2);
    answer.at(3) = shearareaz * g * reducedStrain.at(3);

    return answer;
}

Interface *
LatticeFrameSteelPlastic::giveInterface(InterfaceType type)
{
    return nullptr;
}

void
LatticeFrameSteelPlastic::performRegularReturn(FloatArrayF< 4 > &stress,
                                               double yieldValue,
                                               GaussPoint *gp,
                                               TimeStep *tStep) const
{
    //Use material specific status
    auto status = static_cast< LatticeFrameSteelPlasticStatus * >( this->giveStatus(gp) );

    double deltaLambda = 0.;

    auto trialStress = stress;
    auto tempStress = trialStress;

    //initialise unknowns
    FloatArrayF< 5 >unknowns;
    unknowns.at(1) = trialStress.at(1);
    unknowns.at(2) = trialStress.at(2);
    unknowns.at(3) = trialStress.at(3);
    unknowns.at(4) = trialStress.at(4);
    unknowns.at(5) = 0.;

    yieldValue = computeYieldValue(tempStress, gp, tStep);

    //initiate residuals
    FloatArrayF< 5 >residuals;
    residuals.at(5) = yieldValue;
    double normOfResiduals  = 1.; //just to get into the loop
    int iterationCount = 0;
    while ( normOfResiduals > yieldTol ) {
        iterationCount++;
        if ( iterationCount == newtonIter ) {
            status->letTempReturnResultBe(LatticeFrameSteelPlasticStatus::RR_NotConverged);
            return;
        }

        FloatArrayF< 5 >residualsNorm;
        residualsNorm.at(1) = residuals.at(1) / this->nx0;
        residualsNorm.at(2) = residuals.at(2) / this->mx0;
        residualsNorm.at(3) = residuals.at(3) / this->my0;
        residualsNorm.at(4) = residuals.at(4) / this->mz0;
        residualsNorm.at(5) = residuals.at(5);

        normOfResiduals = norm(residualsNorm);

        if ( std::isnan(normOfResiduals) ) {
            status->letTempReturnResultBe(LatticeFrameSteelPlasticStatus::RR_NotConverged);
            return;
        }

        if ( normOfResiduals > yieldTol ) {
            auto jacobian = computeJacobian(tempStress, deltaLambda, gp, tStep);

            auto solution = solve_check(jacobian, residuals);
            if ( solution.first ) {
                unknowns -= solution.second;
            } else {
                status->letTempReturnResultBe(LatticeFrameSteelPlastic::RR_NotConverged);
            }

            unknowns.at(5) = max(unknowns.at(5), 0.); //Keep deltaLambda greater than zero!

            /* Update increments final values and DeltaLambda*/
            tempStress.at(1) = unknowns.at(1);
            tempStress.at(2) = unknowns.at(2);
            tempStress.at(3) = unknowns.at(3);
            tempStress.at(4) = unknowns.at(4);
            deltaLambda = unknowns.at(5);

            /* Compute the fVector*/
            auto FVector = computeFVector(tempStress, gp, tStep);
            double g = this->e / ( 2. * ( 1. + this->nu ) );
            const double area = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveArea();
            const double ik = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIk();
            const double iy = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIy();
            const double iz = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIz();

            residuals.at(1) = tempStress.at(1) - trialStress.at(1) + area * this->e * deltaLambda * FVector.at(1);
            residuals.at(2) = tempStress.at(2) - trialStress.at(2) + ik * g * deltaLambda * FVector.at(2);
            residuals.at(3) = tempStress.at(3) - trialStress.at(3) + iy * this->e * deltaLambda * FVector.at(3);
            residuals.at(4) = tempStress.at(4) - trialStress.at(4) + iz * this->e * deltaLambda * FVector.at(4);
            residuals.at(5) = computeYieldValue(tempStress, gp, tStep);
        }
    }

    status->letTempReturnResultBe(LatticeFrameSteelPlastic::RR_Converged);
    ;

    stress = tempStress;
}

FloatMatrixF< 5, 5 >
LatticeFrameSteelPlastic::computeJacobian(const FloatArrayF< 4 > &stress,
                                          const double deltaLambda,
                                          GaussPoint *gp,
                                          TimeStep *tStep) const
{
    auto dMMatrix = computeDMMatrix(stress, gp, tStep);
    auto fVector = computeFVector(stress, gp, tStep);
    double g = this->e / ( 2. * ( 1. + this->nu ) );
    const double area = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveArea();
    const double ik = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIk();
    const double iy = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIy();
    const double iz = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIz();


    /* Compute matrix*/
    FloatMatrixF< 5, 5 >jacobian;
    jacobian.at(1, 1) = 1. + this->e * area * deltaLambda * dMMatrix.at(1, 1);
    jacobian.at(1, 2) = 0.;
    jacobian.at(1, 3) = 0.;
    jacobian.at(1, 4) = 0.;
    jacobian.at(1, 5) = this->e * area * fVector.at(1);

    jacobian.at(2, 1) = 0.;
    jacobian.at(2, 2) = 1. + ik * g * deltaLambda * dMMatrix.at(2, 2);
    jacobian.at(2, 3) = 0.;
    jacobian.at(2, 4) = 0.;
    jacobian.at(2, 5) = ik * this->e * fVector.at(2);

    jacobian.at(3, 1) = 0.;
    jacobian.at(3, 2) = 0.;
    jacobian.at(3, 3) = 1. + iy * this->e * deltaLambda * dMMatrix.at(3, 3);
    jacobian.at(3, 4) = 0.;
    jacobian.at(3, 5) = iy * this->e * fVector.at(3);

    jacobian.at(4, 1) = 0.;
    jacobian.at(4, 2) = 0.;
    jacobian.at(4, 3) = 0.;
    jacobian.at(4, 4) = 1. + iz * this->e * deltaLambda * dMMatrix.at(4, 4);
    jacobian.at(4, 5) = iz * this->e * fVector.at(4);

    jacobian.at(5, 1) = fVector.at(1);
    jacobian.at(5, 2) = fVector.at(2);
    jacobian.at(5, 3) = fVector.at(3);
    jacobian.at(5, 4) = fVector.at(4);
    jacobian.at(5, 5) = 0.;

    return jacobian;
}

FloatArrayF< 6 >
LatticeFrameSteelPlastic::giveFrameForces3d(const FloatArrayF< 6 > &originalStrain, GaussPoint *gp, TimeStep *tStep)
{
    auto status = static_cast< LatticeFrameSteelPlasticStatus * >( this->giveStatus(gp) );
    auto reducedStrain = originalStrain;
    auto thermalStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);
    if ( thermalStrain.giveSize() ) {
        reducedStrain -= FloatArrayF< 6 >(thermalStrain);
    }

    auto stress = this->performPlasticityReturn(gp, reducedStrain, tStep);

    status->letTempLatticeStrainBe(originalStrain);
    status->letTempReducedLatticeStrainBe(reducedStrain);
    status->letTempLatticeStressBe(stress);

    return stress;
}



FloatMatrixF< 6, 6 >
LatticeFrameSteelPlastic::give3dFrameStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime) const
{
    static_cast< LatticeFrameSteelPlasticStatus * >( this->giveStatus(gp) );

    double g = this->e / ( 2. * ( 1. + this->nu ) );

    const double area = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveArea();
    const double ik = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIk();
    const double iy = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIy();
    const double iz = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIz();
    const double shearareay = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveShearAreaY();
    const double shearareaz = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveShearAreaZ();

    FloatArrayF< 6 >d = {
        this->e * area,
        g *shearareay,
        g *shearareaz,
        g *ik,
        this->e * iy,
        this->e * iz
    };

    return diag(d);
}

LatticeFrameSteelPlasticStatus::LatticeFrameSteelPlasticStatus(int n, Domain *d, GaussPoint *g) :  LatticeMaterialStatus(g)
{ }

void
LatticeFrameSteelPlasticStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    LatticeMaterialStatus::printOutputAt(file, tStep);

    fprintf(file, "plasticStrains ");
    for ( double s : this->plasticLatticeStrain ) {
        fprintf(file, "% .8e ", s);
    }
}
}
