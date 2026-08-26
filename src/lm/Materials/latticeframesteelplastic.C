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
 *               Copyright (C) 1993 - 2026   Borek Patzak
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
#include "CrossSections/latticecrosssection.h"
#include "engngm.h"
#include "mathfem.h"
#include "Elements/latticestructuralelement.h"
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
LatticeFrameSteelPlastic::initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    LatticeFrameElastic::initializeFrom(ir);

    // Optional automatic capacity definition from yield stress
    this->sigmay = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->sigmay, _IFT_LatticeFrameSteelPlastic_sigmay);

    this->autoCapacitiesFromSigmay = (this->sigmay > 0.0);

    if (!this->autoCapacitiesFromSigmay) {
        IR_GIVE_FIELD(ir, this->nx0, _IFT_LatticeFrameSteelPlastic_nx0);
        IR_GIVE_FIELD(ir, this->mx0, _IFT_LatticeFrameSteelPlastic_mx0);
        IR_GIVE_FIELD(ir, this->my0, _IFT_LatticeFrameSteelPlastic_my0);
        IR_GIVE_FIELD(ir, this->mz0, _IFT_LatticeFrameSteelPlastic_mz0);
    } else {
        // placeholders; actual values computed per element in giveEffectiveCapacities()
        this->nx0 = 0.0;
        this->mx0 = 0.0;
        this->my0 = 0.0;
        this->mz0 = 0.0;
    }

    yieldTol = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->yieldTol, _IFT_LatticeFrameSteelPlastic_tol);

    this->newtonIter = 100;
    IR_GIVE_OPTIONAL_FIELD(ir, this->newtonIter, _IFT_LatticeFrameSteelPlastic_iter);

    numberOfSubIncrements = 10;
    IR_GIVE_OPTIONAL_FIELD(ir, this->numberOfSubIncrements, _IFT_LatticeFrameSteelPlastic_sub);

    this->hardeningLength = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->hardeningLength, _IFT_LatticeFrameSteelPlastic_hlength);

    hType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->hType, _IFT_LatticeFrameSteelPlastic_htype);

    if ( hType == 0 ) {
        H = 0.;
        IR_GIVE_OPTIONAL_FIELD(ir, H, _IFT_LatticeFrameSteelPlastic_h);
    } else if ( hType == 1 ) {
        IR_GIVE_FIELD(ir, h_eps, _IFT_LatticeFrameSteelPlastic_h_eps);
        IR_GIVE_FIELD(ir, h_function_eps, _IFT_LatticeFrameSteelPlastic_h_function_eps);

        if ( h_eps.at(1) != 0. ) {
            throw ValueInputException(ir, _IFT_LatticeFrameSteelPlastic_h_eps,
                                      "The first entry in h_eps must be 0.");
        }

        if ( h_eps.giveSize() != h_function_eps.giveSize() ) {
            throw ValueInputException(ir, _IFT_LatticeFrameSteelPlastic_h_function_eps,
                                      "the size of 'h_eps' and 'h(eps)' must be the same");
        }
    } else {
        throw ValueInputException(ir, _IFT_LatticeFrameSteelPlastic_htype,
                                  "Unknown htype. Should be either 0 or 1.\n");
    }

    capacityMode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, capacityMode, _IFT_LatticeFrameSteelPlastic_capmode);
    if (capacityMode != 0 && capacityMode != 1) {
        throw ValueInputException(ir, _IFT_LatticeFrameSteelPlastic_capmode,
                                  "capmode must be 0 (resultant, default) or 1 (per-width).");
    }
}



void
LatticeFrameSteelPlastic::giveBaseCapacitiesFromSection(double &nx0Base,
                                                        double &mx0Base,
                                                        double &my0Base,
                                                        double &mz0Base,
                                                        GaussPoint *gp) const
{
    if (!this->autoCapacitiesFromSigmay) {
        nx0Base = this->nx0;
        mx0Base = this->mx0;
        my0Base = this->my0;
        mz0Base = this->mz0;
        return;
    }

    auto *elem = static_cast< LatticeStructuralElement * >(gp->giveElement());
    if (!elem) {
        OOFEM_ERROR("Expected LatticeStructuralElement in giveBaseCapacitiesFromSection.");
    }

    const double A  = elem->giveArea(gp);
    const double I1 = elem->giveI1(gp);
    const double I2 = elem->giveI2(gp);
    const double J = elem->giveJ(gp);

    auto *cs = dynamic_cast< LatticeCrossSection * >(elem->giveCrossSection());
    if (!cs) {
        OOFEM_ERROR("Expected LatticeCrossSection in giveBaseCapacitiesFromSection.");
    }

    // ------------------------------------------------------------
    // Circular section
    // ------------------------------------------------------------
    if (cs->giveShape() == 1) {
        const double r = cs->giveRadius();
        if (r <= 0.0) {
            OOFEM_ERROR("Non-positive radius for circular section.");
        }

        // Full plastic axial capacity
        nx0Base = this->sigmay * A;

        // Torsional capacity (this is only a approximation as before)
        mx0Base = this->sigmay * J / r;

        // Full plastic bending capacity of a solid circle:
        const double mpl = this->sigmay * (4.0 * pow(r, 3.0) / 3.0);
        my0Base = mpl;
        mz0Base = mpl;

        return;
    }

    // ------------------------------------------------------------
    // Rectangular section reconstructed from polygon coordinates
    // ------------------------------------------------------------
    double by = 0.0, bz = 0.0;
    if (elem->giveRectangularSectionDimensions(by, bz, gp)) {
        if (by <= 0.0 || bz <= 0.0) {
            OOFEM_ERROR("Non-positive rectangular dimensions in giveBaseCapacitiesFromSection.");
        }

        const double Arect = by * bz;
        if (fabs(Arect - A) > 1e-6 * std::max(1.0, A)) {
            OOFEM_WARNING("Rectangular dimensions from polygon are slightly inconsistent with computed area.");
        }

        // Optional consistency check with principal second moments
        const double Iyy_rect = by * pow(bz, 3.0) / 12.0;
        const double Izz_rect = bz * pow(by, 3.0) / 12.0;

        if (fabs(Iyy_rect - I1) > 1e-5 * std::max(1.0, fabs(I1)) &&
            fabs(Izz_rect - I2) > 1e-5 * std::max(1.0, fabs(I2))) {
            OOFEM_WARNING("Rectangular dimensions are not fully consistent with principal inertia I1/I2.");
        }

        // Full plastic axial capacity
        nx0Base = this->sigmay * Arect;

        // Torsional capacity (simple approximation)
        const double req = 0.5 * std::max(by, bz);
        mx0Base = this->sigmay * J / req;

        // Full plastic bending capacities about local principal section axes
        my0Base = this->sigmay * by * bz * bz / 4.0;
        mz0Base = this->sigmay * bz * by * by / 4.0;

        return;
    }

    OOFEM_ERROR("Automatic sigmay-based capacities currently implemented only for circular sections or rectangular polygon sections.");
}


 
std::unique_ptr< MaterialStatus >
LatticeFrameSteelPlastic::CreateStatus(GaussPoint *gp) const
{
    return std::make_unique< LatticeFrameSteelPlasticStatus >(1, LatticeFrameSteelPlastic::domain, gp);
}

MaterialStatus *
LatticeFrameSteelPlastic::giveStatus(GaussPoint *gp) const
{
    if ( gp->hasMaterialStatus() ) {
        return static_cast < MaterialStatus * > ( gp->giveMaterialStatus() );
    }
    // create a new one
    return static_cast < MaterialStatus * > ( gp->setMaterialStatus(this->CreateStatus(gp)) );
}




double
LatticeFrameSteelPlastic::computeYieldValue(const FloatArrayF < 4 > & stress,
                                            const double kappa,
                                            GaussPoint *gp,
                                            TimeStep *tStep) const
{
    double yieldValue = 0.;
    double hardening = computeHardening(kappa);

    double nx = stress.at(1);
    double mx = stress.at(2);
    double my = stress.at(3);
    double mz = stress.at(4);

    double nx0Eff, mx0Eff, my0Eff, mz0Eff;
    giveEffectiveCapacities(nx0Eff, mx0Eff, my0Eff, mz0Eff, gp, tStep);

    yieldValue = 1. / ( pow(hardening, 2.) ) * ( pow(nx / nx0Eff, 2.) + pow(mx / mx0Eff, 2.) + pow(my / my0Eff, 2.) + pow(mz / mz0Eff, 2.) ) - 1.;

    return yieldValue;
}

void
LatticeFrameSteelPlastic::giveEffectiveCapacities(double &nx0Eff,
                                                  double &mx0Eff,
                                                  double &my0Eff,
                                                  double &mz0Eff,
                                                  GaussPoint *gp,
                                                  TimeStep *tStep) const
{
    double reductionFactor = 1.0;
    if (this->tCrit != 0.0) {
        reductionFactor = computeTemperatureReductionFactor(gp, tStep, VM_Total);
    }

    double widthFactor = 1.0;
    if (this->capacityMode == 1) {
        widthFactor = static_cast< LatticeStructuralElement * >(gp->giveElement())->giveTributaryWidth(gp);
        if (widthFactor <= 0.0) {
            OOFEM_ERROR("Non-positive tributary width.");
        }
    }

    double nx0Base, mx0Base, my0Base, mz0Base;
    giveBaseCapacitiesFromSection(nx0Base, mx0Base, my0Base, mz0Base, gp);

    nx0Eff = reductionFactor * widthFactor * nx0Base;
    mx0Eff = reductionFactor * widthFactor * mx0Base;
    my0Eff = reductionFactor * widthFactor * my0Base;
    mz0Eff = reductionFactor * widthFactor * mz0Base;
}

    
double
LatticeFrameSteelPlastic::computeHardening(const double kappa) const
{
    double hardening = 1.;

    if ( this->hType == 0 ) {
        hardening = 1. + this->H * kappa;
    } else {
        if ( kappa > h_eps.at(h_eps.giveSize() ) ) {
            OOFEM_ERROR("kappa=%e outside range of specified hardening law/n",kappa);
        }

        for ( int i = 1; i < h_eps.giveSize(); i++ ) {
            if ( kappa >= h_eps.at(i) && kappa < h_eps.at(i + 1) ) {
                hardening = h_function_eps.at(i) + ( kappa - h_eps.at(i) ) / ( h_eps.at(i + 1) - h_eps.at(i) ) * ( h_function_eps.at(i + 1) - h_function_eps.at(i) );
                return hardening;
            }
        }
    }

    return hardening;
}


FloatArrayF < 5 >
LatticeFrameSteelPlastic::computeFVector(const FloatArrayF < 4 > & stress,
                                         const double kappa,
                                         GaussPoint * gp,
                                         TimeStep * tStep) const
{
    double hardening = computeHardening(kappa);
    double dHardeningDKappa = computeDHardeningDKappa(kappa);

    double nx0Eff, mx0Eff, my0Eff, mz0Eff;
    giveEffectiveCapacities(nx0Eff, mx0Eff, my0Eff, mz0Eff, gp, tStep);

    double nx = stress.at(1);
    double mx = stress.at(2);
    double my = stress.at(3);
    double mz = stress.at(4);

    FloatArrayF < 5 > f;

    f.at(1) = 2. * nx / ( pow(nx0Eff, 2.) * pow(hardening, 2.) );
    f.at(2) = 2. * mx / ( pow(mx0Eff, 2.) * pow(hardening, 2.) );
    f.at(3) = 2. * my / ( pow(my0Eff, 2.) * pow(hardening, 2.) );
    f.at(4) = 2. * mz / ( pow(mz0Eff, 2.) * pow(hardening, 2.) );
    f.at(5) = -2 / pow(hardening, 3.) * ( pow(nx / nx0Eff, 2.) + pow(mx / mx0Eff, 2.) + pow(my / my0Eff, 2.) + pow(mz / mz0Eff, 2.) ) * dHardeningDKappa;

    return f;
}


FloatArrayF < 5 >
LatticeFrameSteelPlastic::computeMVector(const FloatArrayF < 4 > & stress,
                                         const double kappa,
                                         GaussPoint * gp,
                                         TimeStep * tStep) const
{
    double hardening = computeHardening(kappa);

    double nx0Eff, mx0Eff, my0Eff, mz0Eff;
    giveEffectiveCapacities(nx0Eff, mx0Eff, my0Eff, mz0Eff, gp, tStep);

    double nx = stress.at(1);
    double mx = stress.at(2);
    double my = stress.at(3);
    double mz = stress.at(4);

    FloatArrayF < 5 > m;

    m.at(1) = 2. * nx / ( pow(nx0Eff, 2.) * pow(hardening, 2.) );
    m.at(2) = 2. * mx / ( pow(mx0Eff, 2.) * pow(hardening, 2.) );
    m.at(3) = 2. * my / ( pow(my0Eff, 2.) * pow(hardening, 2.) );
    m.at(4) = 2. * mz / ( pow(mz0Eff, 2.) * pow(hardening, 2.) );

    m.at(5) = sqrt(pow(m.at(1), 2.) + pow(this->hardeningLength, 2.) * ( pow(m.at(2), 2.) + pow(m.at(3), 2.) +  pow(m.at(4), 2.) ) );

    return m;
}




FloatMatrixF < 5, 5 >
LatticeFrameSteelPlastic::computeDMMatrix(const FloatArrayF < 4 > & stress, const double kappa, GaussPoint * gp, TimeStep * tStep) const
{
    auto mVector = computeMVector(stress, kappa, gp, tStep);

    FloatMatrixF < 5, 5 > dm;

    double hardening = computeHardening(kappa);
    double dHardeningDKappa = computeDHardeningDKappa(kappa);

    double nx0Eff, mx0Eff, my0Eff, mz0Eff;
    giveEffectiveCapacities(nx0Eff, mx0Eff, my0Eff, mz0Eff, gp, tStep);

    //dm1dnx
    dm.at(1, 1) = 2. / ( pow(nx0Eff, 2.) * pow(hardening, 2.) );
    //dm2dnx
    dm.at(1, 2) = 0;
    //dm3dnx
    dm.at(1, 3) = 0;
    //dm4dnx
    dm.at(1, 4) = 0;
    //dm5dnx
    dm.at(1, 5) = -4. * stress.at(1) / ( pow(nx0Eff, 2.) * pow(hardening, 3.) ) * dHardeningDKappa;


    dm.at(2, 1) = 0;
    dm.at(2, 2) = 2. / ( pow(mx0Eff, 2.) * pow(hardening, 2.) );
    dm.at(2, 3) = 0;
    dm.at(2, 4) = 0;
    dm.at(2, 5) = -4. * stress.at(2) / ( pow(mx0Eff, 2.) * pow(hardening, 3.) ) * dHardeningDKappa;

    dm.at(3, 1) = 0;
    dm.at(3, 2) = 0;
    dm.at(3, 3) = 2. / ( pow(my0Eff, 2.) * pow(hardening, 2.) );
    dm.at(3, 4) = 0;
    dm.at(3, 5) = -4. * stress.at(3) / ( pow(my0Eff, 2.) * pow(hardening, 3.) ) * dHardeningDKappa;


    dm.at(4, 1) = 0;
    dm.at(4, 2) = 0;
    dm.at(4, 3) = 0;
    dm.at(4, 4) = 2. / ( pow(mz0Eff, 2.) * pow(hardening, 2.) );
    dm.at(4, 5) = -4. * stress.at(4) / ( pow(mz0Eff, 2.) * pow(hardening, 3.) ) * dHardeningDKappa;


    if ( this->hardeningLength == 0. ) {
        dm.at(5, 1) = sgn(stress.at(1) ) * dm.at(1, 1);
        dm.at(5, 2) = 0.0;
        dm.at(5, 3) = 0.0;
        dm.at(5, 4) = 0.0;
        dm.at(5, 5) = -4.0 * fabs(stress.at(1) ) / ( pow(nx0Eff, 2.) * pow(hardening, 3.) ) * dHardeningDKappa;
    } else if ( fabs(mVector.at(5) ) > 1.e-14 )         {
        dm.at(5, 1) = mVector.at(1) * dm.at(1, 1) / mVector.at(5);
        dm.at(5, 2) = mVector.at(2) * dm.at(2, 2) * pow(hardeningLength, 2.) / mVector.at(5);
        dm.at(5, 3) = mVector.at(3) * dm.at(3, 3) * pow(hardeningLength, 2.) / mVector.at(5);
        dm.at(5, 4) = mVector.at(4) * dm.at(4, 4) * pow(hardeningLength, 2.) / mVector.at(5);

        double term1 = mVector.at(1) * ( -4. * stress.at(1) / ( pow(nx0Eff, 2.) * pow(hardening, 3.) ) ) * dHardeningDKappa;
        double term2 = mVector.at(2) * ( -4. * stress.at(2) / ( pow(mx0Eff, 2.) * pow(hardening, 3.) ) ) * dHardeningDKappa;
        double term3 = mVector.at(3) * ( -4. * stress.at(3) / ( pow(my0Eff, 2.) * pow(hardening, 3.) ) ) * dHardeningDKappa;
        double term4 = mVector.at(4) * ( -4. * stress.at(4) / ( pow(mz0Eff, 2.) * pow(hardening, 3.) ) ) * dHardeningDKappa;

        dm.at(5, 5) = ( term1 + pow(this->hardeningLength, 2.) * ( term2 + term3 + term4 ) ) / mVector.at(5);
    } else   {
        OOFEM_ERROR("mVector.at(5) in computeDMMatrix should not be zero\n");
    }
    return dm;
}

FloatArrayF < 6 >
LatticeFrameSteelPlastic::giveThermalDilatationVector(GaussPoint * gp,  TimeStep * tStep) const
// returns a FloatArray(6) of initial strain vector caused by unit temperature in direction of gp (element) local axes
{
    double alpha = this->give(tAlpha, gp);

    return {
        alpha, 0., 0., 0., 0., 0.
    };
}

FloatArrayF < 6 >
LatticeFrameSteelPlastic::giveReducedStrain(GaussPoint * gp, TimeStep * tStep) const
{
    auto status = static_cast < LatticeMaterialStatus * > ( this->giveStatus(gp) );
    return status->giveReducedLatticeStrain();
}


FloatArrayF < 6 >
LatticeFrameSteelPlastic::performPlasticityReturn(GaussPoint * gp, const FloatArrayF < 6 > & reducedStrain, TimeStep * tStep) const
{
    auto status = static_cast < LatticeFrameSteelPlasticStatus * > ( this->giveStatus(gp) );

    //Shear components are not used for plasticity return
    auto strain = reducedStrain [ { 0, 3, 4, 5 } ];

    /* Get plastic strain vector from status*/
    auto tempPlasticStrain = status->givePlasticLatticeStrain() [ { 0, 3, 4, 5 } ];
    auto tempKappa = status->giveKappa();

    FloatArrayF < 6 > tangentFull = diag(give3dFrameStiffnessMatrix(ElasticStiffness, gp, tStep) );

    FloatArrayF < 4 > tangent = {
        tangentFull.at(1), tangentFull.at(4), tangentFull.at(5), tangentFull.at(6)
    };

    /* Compute trial stress*/
    auto stress = mult(tangent, strain - tempPlasticStrain);

    auto oldStrain = this->giveReducedStrain(gp, tStep) [ { 0, 3, 4, 5 } ];

    /* Compute yield value*/
    double yieldValue = computeYieldValue(stress, tempKappa, gp, tStep);
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
            tempKappa = performRegularReturn(stress, yieldValue, gp, tStep);

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
                tempPlasticStrain.at(1) = tempStrain.at(1) - stress.at(1) / ( tangent.at(1) );
                tempPlasticStrain.at(2) = tempStrain.at(2) - stress.at(2) / ( tangent.at(2) );
                tempPlasticStrain.at(3) = tempStrain.at(3) - stress.at(3) / ( tangent.at(3) );
                tempPlasticStrain.at(4) = tempStrain.at(4) - stress.at(4) / ( tangent.at(4) );

                status->letTempPlasticLatticeStrainBe(assemble < 6 > ( tempPlasticStrain, { 0, 3, 4, 5 } ) );

                status->setTempKappa(tempKappa);

                subIncrementFlag = 0;

                status->letTempReturnResultBe(LatticeFrameSteelPlastic::RR_NotConverged);
                convergedStrain = tempStrain;
                deltaStrain = strain - convergedStrain;
                tempStrain = strain;
                subIncrementCounter = 0;
            }
        }
    }


    tempPlasticStrain.at(1) = strain.at(1) - stress.at(1) / ( tangent.at(1) );
    tempPlasticStrain.at(2) = strain.at(2) - stress.at(2) / ( tangent.at(2) );
    tempPlasticStrain.at(3) = strain.at(3) - stress.at(3) / ( tangent.at(3) );
    tempPlasticStrain.at(4) = strain.at(4) - stress.at(4) / ( tangent.at(4) );


    status->letTempPlasticLatticeStrainBe(assemble < 6 > ( tempPlasticStrain, { 0, 3, 4, 5 } ) );

    status->setTempKappa(tempKappa);

    auto answer = assemble < 6 > ( stress, { 0, 3, 4, 5 } );
    answer.at(2) = tangentFull.at(2) * reducedStrain.at(2);
    answer.at(3) = tangentFull.at(3) * reducedStrain.at(3);

    return answer;
}

Interface *
LatticeFrameSteelPlastic::giveInterface(InterfaceType type)
{
    return nullptr;
}

double LatticeFrameSteelPlastic::computeDHardeningDKappa(const double kappa) const
{
    double dHardeningDKappa = 0.;
    if ( hType == 0 ) {
        dHardeningDKappa = this->H;
    } else {
        if ( kappa > h_eps.at(h_eps.giveSize() ) ) {
            OOFEM_ERROR("kappa outside range of specified hardening law/n");
        }
        for ( int i = 1; i < h_eps.giveSize(); i++ ) {
            if ( kappa >= h_eps.at(i) && kappa < h_eps.at(i + 1) ) {
                dHardeningDKappa = ( h_function_eps.at(i + 1) - h_function_eps.at(i) ) / ( h_eps.at(i + 1) - h_eps.at(i) );
            }
        }
    }

    return dHardeningDKappa;
}

double
LatticeFrameSteelPlastic::performRegularReturn(FloatArrayF < 4 > & stress,
                                               double yieldValue,
                                               GaussPoint *gp,
                                               TimeStep *tStep) const
{
    //Use material specific status
    auto status = static_cast < LatticeFrameSteelPlasticStatus * > ( this->giveStatus(gp) );

    FloatArrayF < 6 > tangentFull = diag(give3dFrameStiffnessMatrix(ElasticStiffness, gp, tStep) );

    // Do the same for kappa
    double kappa = status->giveKappa();
    double tempKappa = kappa;


    double deltaLambda = 0.;

    double nx0Eff, mx0Eff, my0Eff, mz0Eff;
    giveEffectiveCapacities(nx0Eff, mx0Eff, my0Eff, mz0Eff, gp, tStep);

    auto trialStress = stress;
    auto tempStress = trialStress;

    //initialise unknowns
    FloatArrayF < 6 > unknowns;
    unknowns.at(1) = trialStress.at(1);
    unknowns.at(2) = trialStress.at(2);
    unknowns.at(3) = trialStress.at(3);
    unknowns.at(4) = trialStress.at(4);
    unknowns.at(5) = tempKappa;
    unknowns.at(6) = 0.;

    yieldValue = computeYieldValue(tempStress, tempKappa, gp, tStep);

    //initiate residuals
    FloatArrayF < 6 > residuals;
    residuals.at(6) = yieldValue;
    double normOfResiduals  = 1.;//just to get into the loop
    int iterationCount = 0;
    while ( normOfResiduals > yieldTol ) {
        iterationCount++;
        if ( iterationCount == newtonIter ) {
            status->letTempReturnResultBe(LatticeFrameSteelPlasticStatus::RR_NotConverged);
            return 0.;
        }

        FloatArrayF < 6 > residualsNorm;
        residualsNorm.at(1) = residuals.at(1) / nx0Eff;
        residualsNorm.at(2) = residuals.at(2) / mx0Eff;
        residualsNorm.at(3) = residuals.at(3) / my0Eff;
        residualsNorm.at(4) = residuals.at(4) / mz0Eff;
        residualsNorm.at(5) = residuals.at(5);
        residualsNorm.at(6) = residuals.at(6);

        normOfResiduals = norm(residualsNorm);

        if ( std::isnan(normOfResiduals) ) {
            status->letTempReturnResultBe(LatticeFrameSteelPlasticStatus::RR_NotConverged);
            return 0.;
        }

        if ( normOfResiduals > yieldTol ) {
            auto jacobian = computeJacobian(tempStress, deltaLambda, tempKappa, gp, tStep);

            auto solution = solve_check(jacobian, residuals);
            if ( solution.first ) {
                unknowns -= solution.second;
            } else {
                status->letTempReturnResultBe(LatticeFrameSteelPlastic::RR_NotConverged);
            }

            unknowns.at(5) = max(unknowns.at(5), kappa);       // Keep deltaKappa greater than zero!
            unknowns.at(6) = max(unknowns.at(6), 0.);       // Keep deltaLambda greater than zero!

            /* Update increments final values and DeltaLambda*/
            tempStress.at(1) = unknowns.at(1);
            tempStress.at(2) = unknowns.at(2);
            tempStress.at(3) = unknowns.at(3);
            tempStress.at(4) = unknowns.at(4);
            tempKappa          = unknowns.at(5);
            deltaLambda        = unknowns.at(6);

            /* Compute the fVector*/
            auto mVector = computeMVector(tempStress, tempKappa, gp, tStep);

            residuals.at(1) = tempStress.at(1) - trialStress.at(1) + tangentFull.at(1) * deltaLambda * mVector.at(1);
            residuals.at(2) = tempStress.at(2) - trialStress.at(2) + tangentFull.at(4) * deltaLambda * mVector.at(2);
            residuals.at(3) = tempStress.at(3) - trialStress.at(3) + tangentFull.at(5) * deltaLambda * mVector.at(3);
            residuals.at(4) = tempStress.at(4) - trialStress.at(4) + tangentFull.at(6) * deltaLambda * mVector.at(4);
            residuals.at(5) = -tempKappa + status->giveKappa() + deltaLambda * mVector.at(5);
            residuals.at(6) = computeYieldValue(tempStress, tempKappa, gp, tStep);
        }
    }
    status->letTempReturnResultBe(LatticeFrameSteelPlastic::RR_Converged);

    stress = tempStress;
    return tempKappa;
}



FloatMatrixF < 6, 6 >
LatticeFrameSteelPlastic::computeJacobian(const FloatArrayF < 4 > & stress,
                                          const double deltaLambda,
                                          const double kappa,
                                          GaussPoint * gp,
                                          TimeStep * tStep) const
{
    auto dMMatrix = computeDMMatrix(stress, kappa, gp, tStep);
    auto fVector = computeFVector(stress, kappa, gp, tStep);
    auto mVector = computeMVector(stress, kappa, gp, tStep);

    FloatArrayF < 6 > tangentFull = diag(give3dFrameStiffnessMatrix(ElasticStiffness, gp, tStep) );

    /* Compute matrix*/
    FloatMatrixF < 6, 6 > jacobian;
    jacobian.at(1, 1) = 1. + tangentFull.at(1) * deltaLambda * dMMatrix.at(1, 1);
    jacobian.at(1, 2) = 0.;
    jacobian.at(1, 3) = 0.;
    jacobian.at(1, 4) = 0.;
    jacobian.at(1, 5) = deltaLambda * tangentFull.at(1) * dMMatrix.at(1, 5);
    jacobian.at(1, 6) = tangentFull.at(1) * mVector.at(1);

    jacobian.at(2, 1) = 0.;
    jacobian.at(2, 2) = 1. + tangentFull.at(4) * deltaLambda * dMMatrix.at(2, 2);
    jacobian.at(2, 3) = 0.;
    jacobian.at(2, 4) = 0.;
    jacobian.at(2, 5) = deltaLambda * tangentFull.at(4) * dMMatrix.at(2, 5);
    jacobian.at(2, 6) = tangentFull.at(4) * mVector.at(2);

    jacobian.at(3, 1) = 0.;
    jacobian.at(3, 2) = 0.;
    jacobian.at(3, 3) = 1. + tangentFull.at(5) * deltaLambda * dMMatrix.at(3, 3);
    jacobian.at(3, 4) = 0.;
    jacobian.at(3, 5) =  deltaLambda * tangentFull.at(5) * dMMatrix.at(3, 5);
    ;
    jacobian.at(3, 6) = tangentFull.at(5) * mVector.at(3);

    jacobian.at(4, 1) = 0.;
    jacobian.at(4, 2) = 0.;
    jacobian.at(4, 3) = 0.;
    jacobian.at(4, 4) = 1. + tangentFull.at(6) * deltaLambda * dMMatrix.at(4, 4);
    jacobian.at(4, 5) = jacobian.at(4, 5) = deltaLambda * tangentFull.at(6) * dMMatrix.at(4, 5);
    jacobian.at(4, 6) = tangentFull.at(6) * mVector.at(4);

    jacobian.at(5, 1) = deltaLambda * dMMatrix.at(5, 1);
    jacobian.at(5, 2) = deltaLambda * dMMatrix.at(5, 2);
    jacobian.at(5, 3) = deltaLambda * dMMatrix.at(5, 3);
    jacobian.at(5, 4) = deltaLambda * dMMatrix.at(5, 4);
    jacobian.at(5, 5) = deltaLambda * dMMatrix.at(5, 5) - 1.;
    jacobian.at(5, 6) = mVector.at(5);

    jacobian.at(6, 1) = fVector.at(1);
    jacobian.at(6, 2) = fVector.at(2);
    jacobian.at(6, 3) = fVector.at(3);
    jacobian.at(6, 4) = fVector.at(4);
    jacobian.at(6, 5) = fVector.at(5);
    jacobian.at(6, 6) = 0.;

    return jacobian;
}

FloatArrayF < 6 >
LatticeFrameSteelPlastic::giveFrameForces3d(const FloatArrayF < 6 > & originalStrain, GaussPoint * gp, TimeStep * tStep)
{
    auto status = static_cast < LatticeFrameSteelPlasticStatus * > ( this->giveStatus(gp) );

    this->initTempStatus(gp);

    auto reducedStrain = originalStrain;
    auto thermalStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);
    if ( thermalStrain.giveSize() ) {
        reducedStrain -= FloatArrayF < 6 > ( thermalStrain );
    }

    auto stress = this->performPlasticityReturn(gp, reducedStrain, tStep);

    status->letTempLatticeStrainBe(originalStrain);
    status->letTempReducedLatticeStrainBe(reducedStrain);
    status->letTempLatticeStressBe(stress);

    return stress;
}


LatticeFrameSteelPlasticStatus::LatticeFrameSteelPlasticStatus(int n, Domain *d, GaussPoint *g) :  LatticeMaterialStatus(g)
{ }

void
LatticeFrameSteelPlasticStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    LatticeMaterialStatus::printOutputAt(file, tStep);

    fprintf(file, "plasticstrain ");
    for ( double s : this->plasticLatticeStrain ) {
        fprintf(file, "% .8e ", s);
    }
    fprintf(file, "kappa  % .8e ", this->kappa);
    fprintf(file, "\n");
}


void
LatticeFrameSteelPlasticStatus::updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    LatticeMaterialStatus::updateYourself(atTime);
    this->kappa = this->tempKappa;
}

void
LatticeFrameSteelPlasticStatus::initTempStatus()
{
    LatticeMaterialStatus::initTempStatus();
    this->tempKappa = this->kappa;
}
}
