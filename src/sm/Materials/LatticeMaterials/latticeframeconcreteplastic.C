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

#include "latticeframeconcreteplastic.h"
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
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "datastream.h"
#include "staggeredproblem.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(LatticeFrameConcretePlastic);

bool LatticeFrameConcretePlastic::hasMaterialModeCapability(MaterialMode mode) const
{
    return ( mode == _3dLattice );
}


void LatticeFrameConcretePlastic::initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    LatticeStructuralMaterial::initializeFrom(ir);

    // Young's modulus of the material that the beam element is made of
    IR_GIVE_FIELD(ir, this->e, _IFT_LatticeFrameConcretePlastic_e); // Macro

    // Poisson's ratio of the material that the beam element is made of
    IR_GIVE_FIELD(ir, this->nu, _IFT_LatticeFrameConcretePlastic_n); // Macro

    // Optional automatic capacity definition from yield stress
    this->sigmay = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->sigmay, _IFT_LatticeFrameConcretePlastic_sigmay);

    this->autoCapacitiesFromSigmay = (this->sigmay > 0.0);


    if (!this->autoCapacitiesFromSigmay) {

        // Nx0
        IR_GIVE_FIELD(ir, this->nx0, _IFT_LatticeFrameConcretePlastic_nx0); // Macro

        // Nx01
        this->nx01 = nx0;
        IR_GIVE_OPTIONAL_FIELD(ir, this->nx01, _IFT_LatticeFrameConcretePlastic_nx01); // Macro

        // vy0
        IR_GIVE_FIELD(ir, this->vy0, _IFT_LatticeFrameConcretePlastic_vy0); // Macro

        // vy01
        this->vy01 = vy0;
        IR_GIVE_OPTIONAL_FIELD(ir, this->vy01, _IFT_LatticeFrameConcretePlastic_vy01); // Macro

        // vz0
        IR_GIVE_FIELD(ir, this->vz0, _IFT_LatticeFrameConcretePlastic_vz0); // Macro

        // vz01
        this->vz01 = vz0;
        IR_GIVE_OPTIONAL_FIELD(ir, this->vz01, _IFT_LatticeFrameConcretePlastic_vz01); // Macro

        // Mx0
        IR_GIVE_FIELD(ir, this->mx0, _IFT_LatticeFrameConcretePlastic_mx0); // Macro

        // Mx01
        this->mx01 = mx0;
        IR_GIVE_OPTIONAL_FIELD(ir, this->mx01, _IFT_LatticeFrameConcretePlastic_mx01); // Macro

        // My0
        IR_GIVE_FIELD(ir, this->my0, _IFT_LatticeFrameConcretePlastic_my0); // Macro

        // My01
        this->my01 = my0;
        IR_GIVE_OPTIONAL_FIELD(ir, this->my01, _IFT_LatticeFrameConcretePlastic_my01); // Macro

        // Mz0
        IR_GIVE_FIELD(ir, this->mz0, _IFT_LatticeFrameConcretePlastic_mz0); // Macro

        // Mz01
        this->mz01 = mz0;
        IR_GIVE_OPTIONAL_FIELD(ir, this->mz01, _IFT_LatticeFrameConcretePlastic_mz01); // Macro

    }
    else {
        // placeholders; actual values computed per element in giveEffectiveCapacities()
        this->nx0 = 0.0;
        this->vy0 = 0.0;
        this->vz0 = 0.0;
        this->mx0 = 0.0;
        this->my0 = 0.0;
        this->mz0 = 0.0;
        this->nx01 = 0.0;
        this->vy01 = 0.0;
        this->vz01 = 0.0;
        this->mx01 = 0.0;
        this->my01 = 0.0;
        this->mz01 = 0.0;
    }

    yieldTol = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->yieldTol, _IFT_LatticeFrameConcretePlastic_tol); // Macro
    
    this->newtonIter = 100;
    IR_GIVE_OPTIONAL_FIELD(ir, this->newtonIter, _IFT_LatticeFrameConcretePlastic_iter); // Macro

    numberOfSubIncrements = 10;
    IR_GIVE_OPTIONAL_FIELD(ir, this->numberOfSubIncrements, _IFT_LatticeFrameConcretePlastic_sub); // Macro

    this->plasticFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, plasticFlag, _IFT_LatticeFrameConcretePlastic_plastic); // Macro

    // wu
    this->wu = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->wu, _IFT_LatticeFrameConcretePlastic_wu);

    // wf
    this->wf = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, wf, _IFT_LatticeFrameConcretePlastic_wf);

    capacityMode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, capacityMode, _IFT_LatticeFrameConcretePlastic_capmode);
    if (capacityMode != 0 && capacityMode != 1) {
        throw ValueInputException(ir, _IFT_LatticeFrameConcretePlastic_capmode,
                                  "capmode must be 0 (resultant, default) or 1 (per-width).");
    }

}



LatticeFrameConcretePlastic::CapacitySet
LatticeFrameConcretePlastic::giveBaseCapacitiesFromSection(GaussPoint *gp) const
{
    LatticeFrameConcretePlastic::CapacitySet cap;

    if (!this->autoCapacitiesFromSigmay) {
        cap.nx0  = this->nx0;
        cap.vy0  = this->vy0;
        cap.vz0  = this->vz0;
        cap.mx0  = this->mx0;
        cap.my0  = this->my0;
        cap.mz0  = this->mz0;

        cap.nx01 = this->nx01;
        cap.vy01 = this->vy01;
        cap.vz01 = this->vz01;
        cap.mx01 = this->mx01;
        cap.my01 = this->my01;
        cap.mz01 = this->mz01;

        return cap;
    }

    auto *elem = static_cast<LatticeStructuralElement *>(gp->giveElement());
    if (!elem) {
        OOFEM_ERROR("Expected LatticeStructuralElement.");
    }

    const double A  = elem->giveArea(gp);
    const double I1 = elem->giveI1(gp);
    const double I2 = elem->giveI2(gp);
    const double J = elem->giveJ(gp);

    auto *cs = dynamic_cast<LatticeCrossSection *>(elem->giveCrossSection());
    if (!cs) {
        OOFEM_ERROR("Expected LatticeCrossSection.");
    }

    // ------------------ circular ------------------
    if (cs->giveShape() == 1) {
        const double r = cs->giveRadius();
        if (r <= 0.0) {
            OOFEM_ERROR("Non-positive radius.");
        }

        cap.nx0 = this->sigmay * A;

        cap.vy0 = this->sigmay * A;
        cap.vz0 = this->sigmay * A;

        cap.mx0 = this->sigmay * J / r;

        const double mpl = this->sigmay * (4.0 * pow(r, 3.0) / 3.0);
        cap.my0 = mpl;
        cap.mz0 = mpl;

        // symmetric
        cap.nx01 = cap.nx0;
        cap.vy01 = cap.vy0;
        cap.vz01 = cap.vz0;
        cap.mx01 = cap.mx0;
        cap.my01 = cap.my0;
        cap.mz01 = cap.mz0;

        return cap;
    }

    // ------------------ rectangle ------------------
    double by = 0.0, bz = 0.0;
    if (elem->giveRectangularSectionDimensions(by, bz, gp)) {

        const double Arect = by * bz;

        cap.nx0 = this->sigmay * Arect;

        cap.vy0 = this->sigmay * Arect;
        cap.vz0 = this->sigmay * Arect;

        const double req = 0.5 * std::max(by, bz);
        cap.mx0 = this->sigmay * J / req;

        cap.my0 = this->sigmay * by * bz * bz / 4.0;
        cap.mz0 = this->sigmay * bz * by * by / 4.0;

        // symmetric
        cap.nx01 = cap.nx0;
        cap.vy01 = cap.vy0;
        cap.vz01 = cap.vz0;
        cap.mx01 = cap.mx0;
        cap.my01 = cap.my0;
        cap.mz01 = cap.mz0;

        return cap;
    }

    OOFEM_ERROR("Unsupported section for automatic capacities.");
}


LatticeFrameConcretePlastic::CapacitySet
LatticeFrameConcretePlastic::giveEffectiveCapacities(GaussPoint *gp,
                                                     TimeStep *tStep) const
{
    LatticeFrameConcretePlastic::CapacitySet base = giveBaseCapacitiesFromSection(gp);
    LatticeFrameConcretePlastic::CapacitySet eff;

    double reductionFactor = 1.0;
    if (this->tCrit != 0.0) {
        reductionFactor = computeTemperatureReductionFactor(gp, tStep, VM_Total);
    }

    double widthFactor = 1.0;
    if (this->capacityMode == 1) {
        widthFactor = static_cast<LatticeStructuralElement *>(gp->giveElement())->giveTributaryWidth(gp);
        if (widthFactor <= 0.0) {
            OOFEM_ERROR("Non-positive tributary width.");
        }
    }

    const double f = reductionFactor * widthFactor;

    eff.nx0  = f * base.nx0;
    eff.vy0  = f * base.vy0;
    eff.vz0  = f * base.vz0;
    eff.mx0  = f * base.mx0;
    eff.my0  = f * base.my0;
    eff.mz0  = f * base.mz0;

    eff.nx01 = f * base.nx01;
    eff.vy01 = f * base.vy01;
    eff.vz01 = f * base.vz01;
    eff.mx01 = f * base.mx01;
    eff.my01 = f * base.my01;
    eff.mz01 = f * base.mz01;

    return eff;
}

 
   
std::unique_ptr< MaterialStatus >
LatticeFrameConcretePlastic::CreateStatus(GaussPoint *gp) const
{
    return std::make_unique< LatticeFrameConcretePlasticStatus >(1, LatticeFrameConcretePlastic::domain, gp);
}

MaterialStatus *
LatticeFrameConcretePlastic::giveStatus(GaussPoint *gp) const
{
    if ( gp->hasMaterialStatus() ) {
        return static_cast < MaterialStatus * > ( gp->giveMaterialStatus() );
    }
    // create a new one
    return static_cast < MaterialStatus * > ( gp->setMaterialStatus(this->CreateStatus(gp)) );
}
int
LatticeFrameConcretePlastic::giveIPValue(FloatArray &answer,
                                         GaussPoint *gp,
                                         InternalStateType type,
                                         TimeStep *atTime)
{
    auto status = static_cast < LatticeFrameConcretePlasticStatus * > ( this->giveStatus(gp) );

    if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveDamage();
        return 1;
    } else {
        return LatticeStructuralMaterial::giveIPValue(answer, gp, type, atTime);
    }
}


double
LatticeFrameConcretePlastic::computeYieldValue(const FloatArrayF<6> &stress,
                                               const FloatArrayF<6> &k,
                                               GaussPoint *gp,
                                               TimeStep *tStep) const
{
    const double nx = stress.at(1);
    const double vy = stress.at(2);
    const double vz = stress.at(3);
    const double mx = stress.at(4);
    const double my = stress.at(5);
    const double mz = stress.at(6);

    LatticeFrameConcretePlastic::CapacitySet cap = this->giveEffectiveCapacities(gp, tStep);

    const double ax = (k.at(1) == 1.0) ? cap.nx0  : cap.nx01;
    const double ay = (k.at(2) == 1.0) ? cap.vy0  : cap.vy01;
    const double az = (k.at(3) == 1.0) ? cap.vz0  : cap.vz01;
    const double bx = (k.at(4) == 1.0) ? cap.mx0  : cap.mx01;
    const double by = (k.at(5) == 1.0) ? cap.my0  : cap.my01;
    const double bz = (k.at(6) == 1.0) ? cap.mz0  : cap.mz01;

    if (ax <= 0.0 || ay <= 0.0 || az <= 0.0 || bx <= 0.0 || by <= 0.0 || bz <= 0.0) {
        OOFEM_ERROR("Non-positive effective capacity in computeYieldValue.");
    }

    const double yieldValue =
        pow(nx / ax, 2.0) +
        pow(vy / ay, 2.0) +
        pow(vz / az, 2.0) +
        pow(mx / bx, 2.0) +
        pow(my / by, 2.0) +
        pow(mz / bz, 2.0) - 1.0;

    return yieldValue;
}


FloatArrayF<6>
LatticeFrameConcretePlastic::computeFVector(const FloatArrayF<6> &stress,
                                            const FloatArrayF<6> &k,
                                            GaussPoint *gp,
                                            TimeStep *tStep) const
{
    const double nx = stress.at(1);
    const double vy = stress.at(2);
    const double vz = stress.at(3);
    const double mx = stress.at(4);
    const double my = stress.at(5);
    const double mz = stress.at(6);

    LatticeFrameConcretePlastic::CapacitySet cap = this->giveEffectiveCapacities(gp, tStep);

    const double ax = (k.at(1) == 1.0) ? cap.nx0  : cap.nx01;
    const double ay = (k.at(2) == 1.0) ? cap.vy0  : cap.vy01;
    const double az = (k.at(3) == 1.0) ? cap.vz0  : cap.vz01;
    const double bx = (k.at(4) == 1.0) ? cap.mx0  : cap.mx01;
    const double by = (k.at(5) == 1.0) ? cap.my0  : cap.my01;
    const double bz = (k.at(6) == 1.0) ? cap.mz0  : cap.mz01;

    if (ax <= 0.0 || ay <= 0.0 || az <= 0.0 ||
        bx <= 0.0 || by <= 0.0 || bz <= 0.0) {
        OOFEM_ERROR("Non-positive effective capacity in computeFVector.");
    }

    FloatArrayF<6> f;

    const double ax2 = ax * ax;
    const double ay2 = ay * ay;
    const double az2 = az * az;
    const double bx2 = bx * bx;
    const double by2 = by * by;
    const double bz2 = bz * bz;

    f.at(1) = 2.0 * nx / ax2;
    f.at(2) = 2.0 * vy / ay2;
    f.at(3) = 2.0 * vz / az2;
    f.at(4) = 2.0 * mx / bx2;
    f.at(5) = 2.0 * my / by2;
    f.at(6) = 2.0 * mz / bz2;

    return f;
}
    

FloatMatrixF<6,6>
LatticeFrameConcretePlastic::computeDMMatrix(const FloatArrayF<6> &stress,
                                             const FloatArrayF<6> &k,
                                             GaussPoint *gp,
                                             TimeStep *tStep) const
{
    LatticeFrameConcretePlastic::CapacitySet cap = this->giveEffectiveCapacities(gp, tStep);

    const double ax = (k.at(1) == 1.0) ? cap.nx0  : cap.nx01;
    const double ay = (k.at(2) == 1.0) ? cap.vy0  : cap.vy01;
    const double az = (k.at(3) == 1.0) ? cap.vz0  : cap.vz01;
    const double bx = (k.at(4) == 1.0) ? cap.mx0  : cap.mx01;
    const double by = (k.at(5) == 1.0) ? cap.my0  : cap.my01;
    const double bz = (k.at(6) == 1.0) ? cap.mz0  : cap.mz01;

    if (ax <= 0.0 || ay <= 0.0 || az <= 0.0 ||
        bx <= 0.0 || by <= 0.0 || bz <= 0.0) {
        OOFEM_ERROR("Non-positive effective capacity in computeDMMatrix.");
    }

    const FloatArrayF<6> dm = {
        2.0 / (ax * ax),
        2.0 / (ay * ay),
        2.0 / (az * az),
        2.0 / (bx * bx),
        2.0 / (by * by),
        2.0 / (bz * bz)
    };

    return diag(dm);
}


    
FloatArrayF < 6 >
LatticeFrameConcretePlastic::giveThermalDilatationVector(GaussPoint * gp, TimeStep * tStep) const
// returns a FloatArray(6) of initial strain vector caused by unit temperature in direction of gp (element) local axes
{
    double alpha = this->give(tAlpha, gp);

    return {
        alpha, 0., 0., 0., 0., 0.
    };
}

FloatArrayF < 6 >
LatticeFrameConcretePlastic::giveStrain(GaussPoint * gp, TimeStep * tStep) const
{
    auto status = static_cast < LatticeMaterialStatus * > ( this->giveStatus(gp) );
    return status->giveLatticeStrain();
}

const FloatArrayF < 6 >
LatticeFrameConcretePlastic::checkTransition(const FloatArrayF < 6 > & stress,  GaussPoint * gp, TimeStep * tStep) const
{
    FloatArrayF < 6 > k;

    IntArray m(6);

    if ( stress.at(1) >= 0 ) {
        m.at(1) = 1;
    } else {
        m.at(1) = -1;
    }

    if ( stress.at(2) >= 0 ) {
        m.at(2) = 1;
    } else {
        m.at(2) = -1;
    }

    if ( stress.at(3) >= 0 ) {
        m.at(3) = 1;
    } else {
        m.at(3) = -1;
    }

    if ( stress.at(4) >= 0 ) {
        m.at(4) = 1;
    } else {
        m.at(4) = -1;
    }

    if ( stress.at(5) >= 0 ) {
        m.at(5) = 1;
    } else {
        m.at(5) = -1;
    }

    if ( stress.at(6) >= 0 ) {
        m.at(6) = 1;
    } else {
        m.at(6) = -1;
    }

    int counter = 1;
    for ( int xi = -1; xi < 2; xi += 2 ) {
        for ( int yi = -1; yi < 2; yi += 2 ) {
            for ( int zi = -1; zi < 2; zi += 2 ) {
                for ( int xj = -1; xj < 2; xj += 2 ) {
                    for ( int yj = -1; yj < 2; yj += 2 ) {
                        for ( int zj = -1; zj < 2; zj += 2 ) {
                            if ( !( zj == 0 && yj == 0 && xj == 0 && zi == 0 && yi == 0 && xi == 0 ) ) {
                                if ( xi ==  m.at(1) && yi ==  m.at(2) && zi ==  m.at(3) && xj ==  m.at(4) && yj ==  m.at(5) && zj ==  m.at(6) ) {
                                    k.at(1) = xi;
                                    k.at(2) = yi;
                                    k.at(3) = zi;
                                    k.at(4) = xj;
                                    k.at(5) = yj;
                                    k.at(6) = zj;
                                }
                                counter++;
                            }
                        }
                    }
                }
            }
        }
    }
    if
        ( k.at(1) != m.at(1) && k.at(2) != m.at(2) && k.at(2) != m.at(2) && k.at(2) != m.at(2) && k.at(2) != m.at(2) && k.at(2) != m.at(2) ) {
        OOFEM_ERROR("This case should not exist");
    }

    return k;
}

//

FloatArrayF < 6 >
LatticeFrameConcretePlastic::performPlasticityReturn(GaussPoint * gp, const FloatArrayF < 6 > & strain, TimeStep * tStep) const
{
    LatticeFrameConcretePlastic_ReturnResult returnResult = RR_Unknown;
    int kIter = 0;
    double g = this->e / ( 2. * ( 1. + this->nu ) );
    const double area       = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveArea(gp);
    const double shearareay = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveShearArea1(gp);
    const double shearareaz = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveShearArea2(gp);
    const double ik         = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveJ(gp);
    const double iy         = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveI1(gp);
    const double iz         = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveI2(gp);

    auto status = static_cast < LatticeFrameConcretePlasticStatus * > ( this->giveStatus(gp) );


    /* Get plastic strain vector from status*/
    auto tempPlasticStrain = status->givePlasticLatticeStrain();

    FloatArrayF < 6 > tangent = {
        area * this->e, g * shearareay, g * shearareaz, ik * g, iy * this->e, iz * this->e
    };

    /* Compute trial stress*/
    auto stress = mult(tangent, strain - tempPlasticStrain);



    auto oldStrain = this->giveStrain(gp, tStep);
    //test
    auto k = checkTransition(stress, gp, tStep);

    /* Compute yield value*/
    double yieldValue = computeYieldValue(stress, k, gp, tStep);
    int subIncrementCounter = 0;

    /* Check yield condition, i.e. if the yield value is less than the yield tolerance. If yield condition is valid. Do perform regular return (closest point return)*/

    if ( yieldValue > yieldTol ) {
        int subIncrementFlag = 0;
        auto convergedStrain = oldStrain;
        auto tempStrain      = strain;
        auto deltaStrain     = strain - oldStrain;
        // To get into the loop
        status->letTempReturnResultBe(LatticeFrameConcretePlastic::RR_NotConverged);
        while ( status->giveTempReturnResult() == RR_NotConverged || subIncrementFlag == 1 ) {
            stress = mult(tangent, tempStrain - tempPlasticStrain);
            performRegularReturn(stress, returnResult, k, yieldValue, gp, tStep);
            if ( status->giveTempReturnResult() == RR_NotConverged ) {
                subIncrementCounter++;
                if ( subIncrementCounter > numberOfSubIncrements ) {
                    OOFEM_LOG_INFO("Unstable element %d \n", gp->giveElement()->giveGlobalNumber() );
                    OOFEM_LOG_INFO("Yield value %e \n", yieldValue);
                    OOFEM_LOG_INFO("ConvergedStrain value %e %e %e %e\n", convergedStrain.at(1), convergedStrain.at(2), convergedStrain.at(3), convergedStrain.at(4), convergedStrain.at(5), convergedStrain.at(6) );
                    OOFEM_LOG_INFO("tempStrain value %e %e %e %e\n", tempStrain.at(1), tempStrain.at(2), tempStrain.at(3), tempStrain.at(4), tempStrain.at(5), tempStrain.at(6) );
                    OOFEM_LOG_INFO("deltaStrain value %e %e %e %e\n", deltaStrain.at(1), deltaStrain.at(2), deltaStrain.at(3), deltaStrain.at(4), deltaStrain.at(5), deltaStrain.at(6) );
                    OOFEM_LOG_INFO("targetstrain value %e %e %e %e\n", strain.at(1), strain.at(2), strain.at(3), strain.at(4), strain.at(5), strain.at(6) );

                    OOFEM_ERROR("LatticeFrameConcretePlastic :: performPlasticityReturn - Could not reach convergence with small deltaStrain, giving up.");
                }
                subIncrementFlag = 1;
                deltaStrain *= 0.5;
                tempStrain = convergedStrain + deltaStrain;
            } else if ( status->giveTempReturnResult() == RR_Converged && subIncrementFlag == 1 ) {
                tempPlasticStrain.at(1) = tempStrain.at(1) - stress.at(1) / ( area * this->e );
                tempPlasticStrain.at(2) = tempStrain.at(2) - stress.at(2) / ( g * shearareay );
                tempPlasticStrain.at(3) = tempStrain.at(3) - stress.at(3) / ( g * shearareaz );
                tempPlasticStrain.at(4) = tempStrain.at(4) - stress.at(4) / ( ik * g );
                tempPlasticStrain.at(5) = tempStrain.at(5) - stress.at(5) / ( iy * this->e );
                tempPlasticStrain.at(6) = tempStrain.at(6) - stress.at(6) / ( iz * this->e );

                status->letTempPlasticLatticeStrainBe(assemble < 6 > ( tempPlasticStrain, { 0, 1, 2, 3, 4, 5 } ) );

                subIncrementFlag = 0;

                status->letTempReturnResultBe(LatticeFrameConcretePlastic::RR_NotConverged);
                convergedStrain     = tempStrain;
                deltaStrain         = strain - convergedStrain;
                tempStrain          = strain;
                subIncrementCounter = 0;
            } else { //Converged
                //Check the surface
                FloatArrayF < 6 > kCheck = checkTransition(stress, gp, tStep);

                if   ( kCheck.at(1) == k.at(1) && kCheck.at(2) == k.at(2) && kCheck.at(3) == k.at(3) && kCheck.at(4) == k.at(4) && kCheck.at(5) == k.at(5) &&  kCheck.at(6) == k.at(6) ) {
                    returnResult = RR_Converged;
                } else {
                    //ended up in a region for which other surface should have been used.
                    //Try with other surface
                    if ( kIter == 1 ) {
                        OOFEM_ERROR("LatticePlasticityDamage :: performPlasticityReturn - Tried both surfaces");
                    }

                    k = kCheck;
                    returnResult = RR_NotConverged;
                    kIter++;
                }
            }
        }
    }

    tempPlasticStrain.at(1) = strain.at(1) - stress.at(1) / ( area * this->e );
    tempPlasticStrain.at(2) = strain.at(2) - stress.at(2) / ( g * shearareay );
    tempPlasticStrain.at(3) = strain.at(3) - stress.at(3) / ( g * shearareaz );
    tempPlasticStrain.at(4) = strain.at(4) - stress.at(4) / ( ik * g );
    tempPlasticStrain.at(5) = strain.at(5) - stress.at(5) / ( iy * this->e );
    tempPlasticStrain.at(6) = strain.at(6) - stress.at(6) / ( iz * this->e );

    status->letTempPlasticLatticeStrainBe(tempPlasticStrain);

    return stress;
}

Interface *
LatticeFrameConcretePlastic::giveInterface(InterfaceType type)
{
    return nullptr;
}


void
LatticeFrameConcretePlastic::performRegularReturn(FloatArrayF<6> &stress,
                                                  LatticeFrameConcretePlastic_ReturnResult,
                                                  const FloatArrayF<6> &k,
                                                  double yieldValue,
                                                  GaussPoint *gp,
                                                  TimeStep *tStep) const
{
    // Use material specific status
    auto status = static_cast<LatticeFrameConcretePlasticStatus *>( this->giveStatus(gp) );

    double deltaLambda = 0.0;

    auto trialStress = stress;
    auto tempStress  = trialStress;

    // initialise unknowns
    FloatArrayF<7> unknowns;
    unknowns.at(1) = trialStress.at(1);
    unknowns.at(2) = trialStress.at(2);
    unknowns.at(3) = trialStress.at(3);
    unknowns.at(4) = trialStress.at(4);
    unknowns.at(5) = trialStress.at(5);
    unknowns.at(6) = trialStress.at(6);
    unknowns.at(7) = 0.0;

    yieldValue = computeYieldValue(tempStress, k, gp, tStep);

    // initiate residuals
    FloatArrayF<7> residuals;
    residuals.at(7) = yieldValue;

    double normOfResiduals = 1.0; // just to get into the loop
    int iterationCount = 0;

    while ( normOfResiduals > yieldTol ) {
        iterationCount++;
        if ( iterationCount == newtonIter ) {
            status->letTempReturnResultBe(LatticeFrameConcretePlasticStatus::RR_NotConverged);
            return;
        }

        LatticeFrameConcretePlastic::CapacitySet cap = this->giveEffectiveCapacities(gp, tStep);

        const double ax = (k.at(1) == 1.0) ? cap.nx0  : cap.nx01;
        const double ay = (k.at(2) == 1.0) ? cap.vy0  : cap.vy01;
        const double az = (k.at(3) == 1.0) ? cap.vz0  : cap.vz01;
        const double bx = (k.at(4) == 1.0) ? cap.mx0  : cap.mx01;
        const double by = (k.at(5) == 1.0) ? cap.my0  : cap.my01;
        const double bz = (k.at(6) == 1.0) ? cap.mz0  : cap.mz01;

        if (ax <= 0.0 || ay <= 0.0 || az <= 0.0 ||
            bx <= 0.0 || by <= 0.0 || bz <= 0.0) {
            OOFEM_ERROR("Non-positive effective capacity in performRegularReturn.");
        }

        FloatArrayF<7> residualsNorm;
        residualsNorm.at(1) = residuals.at(1) / ax;
        residualsNorm.at(2) = residuals.at(2) / ay;
        residualsNorm.at(3) = residuals.at(3) / az;
        residualsNorm.at(4) = residuals.at(4) / bx;
        residualsNorm.at(5) = residuals.at(5) / by;
        residualsNorm.at(6) = residuals.at(6) / bz;
        residualsNorm.at(7) = residuals.at(7);

        normOfResiduals = norm(residualsNorm);

        if ( std::isnan(normOfResiduals) ) {
            status->letTempReturnResultBe(LatticeFrameConcretePlasticStatus::RR_NotConverged);
            return;
        }

        if ( normOfResiduals > yieldTol ) {
            auto jacobian = computeJacobian(tempStress, k, deltaLambda, gp, tStep);

            auto solution = solve_check(jacobian, residuals);
            if ( solution.first ) {
                unknowns -= solution.second;
            } else {
                status->letTempReturnResultBe(LatticeFrameConcretePlastic::RR_NotConverged);
            }

            unknowns.at(7) = max(unknowns.at(7), 0.0); // Keep deltaLambda greater than zero!

            // Update unknowns
            tempStress.at(1) = unknowns.at(1);
            tempStress.at(2) = unknowns.at(2);
            tempStress.at(3) = unknowns.at(3);
            tempStress.at(4) = unknowns.at(4);
            tempStress.at(5) = unknowns.at(5);
            tempStress.at(6) = unknowns.at(6);
            deltaLambda = unknowns.at(7);

            // Compute flow vector
            auto FVector = computeFVector(tempStress, k, gp, tStep);

            const double g = this->e / ( 2.0 * ( 1.0 + this->nu ) );

            auto *elem = static_cast<LatticeStructuralElement *>( gp->giveElement() );
            const double area       = elem->giveArea(gp);
            const double shearareay = elem->giveShearArea1(gp);
            const double shearareaz = elem->giveShearArea2(gp);
            const double ik         = elem->giveJ(gp);
            const double iy         = elem->giveI1(gp);
            const double iz         = elem->giveI2(gp);

            residuals.at(1) = tempStress.at(1) - trialStress.at(1) + area       * this->e * deltaLambda * FVector.at(1);
            residuals.at(2) = tempStress.at(2) - trialStress.at(2) + shearareay * g       * deltaLambda * FVector.at(2);
            residuals.at(3) = tempStress.at(3) - trialStress.at(3) + shearareaz * g       * deltaLambda * FVector.at(3);
            residuals.at(4) = tempStress.at(4) - trialStress.at(4) + ik         * g       * deltaLambda * FVector.at(4);
            residuals.at(5) = tempStress.at(5) - trialStress.at(5) + iy         * this->e * deltaLambda * FVector.at(5);
            residuals.at(6) = tempStress.at(6) - trialStress.at(6) + iz         * this->e * deltaLambda * FVector.at(6);
            residuals.at(7) = computeYieldValue(tempStress, k, gp, tStep);
        }
    }

    status->letTempReturnResultBe(LatticeFrameConcretePlastic::RR_Converged);
    stress = tempStress;
}





    
FloatMatrixF < 7, 7 >
LatticeFrameConcretePlastic::computeJacobian(const FloatArrayF < 6 > & stress, const FloatArrayF < 6 > & k,
                                             const double deltaLambda,
                                             GaussPoint * gp,
                                             TimeStep * tStep) const
{
    auto dMMatrix           = computeDMMatrix(stress, k, gp, tStep);
    auto fVector            = computeFVector(stress, k, gp, tStep);
    double g = this->e / ( 2. * ( 1. + this->nu ) );
    const double area       = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveArea(gp);
    const double shearareay = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveShearArea1(gp);
    const double shearareaz = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveShearArea2(gp);
    const double ik         = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveJ(gp);
    const double iy         = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveI1(gp);
    const double iz         = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveI2(gp);

    /* Compute matrix*/
    FloatMatrixF < 7, 7 > jacobian;
    jacobian.at(1, 1) = 1. + this->e * area * deltaLambda * dMMatrix.at(1, 1);
    jacobian.at(1, 2) = 0.;
    jacobian.at(1, 3) = 0.;
    jacobian.at(1, 4) = 0.;
    jacobian.at(1, 5) = 0.;
    jacobian.at(1, 6) = 0.;
    jacobian.at(1, 7) = this->e * area * fVector.at(1);

    jacobian.at(2, 1) = 0.;
    jacobian.at(2, 2) = 1. + shearareay * g * deltaLambda * dMMatrix.at(2, 2);
    jacobian.at(2, 3) = 0.;
    jacobian.at(2, 4) = 0.;
    jacobian.at(2, 5) = 0.;
    jacobian.at(2, 6) = 0.;
    jacobian.at(2, 7) = shearareay * g * fVector.at(2);

    jacobian.at(3, 1) = 0.;
    jacobian.at(3, 2) = 0.;
    jacobian.at(3, 3) = 1. + shearareaz * g * deltaLambda * dMMatrix.at(3, 3);
    jacobian.at(3, 4) = 0.;
    jacobian.at(3, 5) = 0.;
    jacobian.at(3, 6) = 0.;
    jacobian.at(3, 7) = shearareaz * g * fVector.at(3);

    jacobian.at(4, 1) = 0.;
    jacobian.at(4, 2) = 0.;
    jacobian.at(4, 3) = 0.;
    jacobian.at(4, 4) = 1. + ik * g * deltaLambda * dMMatrix.at(4, 4);
    jacobian.at(4, 5) = 0.;
    jacobian.at(4, 6) = 0.;
    jacobian.at(4, 7) = ik * g * fVector.at(4);

    jacobian.at(5, 1) = 0.;
    jacobian.at(5, 2) = 0.;
    jacobian.at(5, 3) = 0.;
    jacobian.at(5, 4) = 0.;
    jacobian.at(5, 5) = 1. + iy * this->e * deltaLambda * dMMatrix.at(5, 5);
    jacobian.at(5, 6) = 0.;
    jacobian.at(5, 7) = iy * this->e * fVector.at(5);


    jacobian.at(6, 1) = 0.;
    jacobian.at(6, 2) = 0.;
    jacobian.at(6, 3) = 0.;
    jacobian.at(6, 4) = 0.;
    jacobian.at(6, 5) = 0.;
    jacobian.at(6, 6) = 1. + this->e * iz * deltaLambda * dMMatrix.at(6, 6);
    jacobian.at(6, 7) = this->e * iz * fVector.at(6);

    jacobian.at(7, 1) = fVector.at(1);
    jacobian.at(7, 2) = fVector.at(2);
    jacobian.at(7, 3) = fVector.at(3);
    jacobian.at(7, 4) = fVector.at(4);
    jacobian.at(7, 5) = fVector.at(5);
    jacobian.at(7, 6) = fVector.at(6);
    jacobian.at(7, 7) = 0.;

    return jacobian;
}


FloatArrayF < 6 >
LatticeFrameConcretePlastic::giveFrameForces3d(const FloatArrayF < 6 > & originalStrain, GaussPoint * gp, TimeStep * tStep)
{
      
      
    auto status        = static_cast < LatticeFrameConcretePlasticStatus * > ( this->giveStatus(gp) );

    this->initTempStatus(gp);
      
    auto strain = originalStrain;
    auto thermalStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);
    if ( thermalStrain.giveSize() ) {
        strain -= FloatArrayF < 6 > ( thermalStrain );
    }
    auto stress = this->performPlasticityReturn(gp, strain, tStep);


    double tempKappaD = 0.0;
    double omega = 0.0;

    //Need to switch off damage if plasticFlag =1
    if(plasticFlag == 0){
        auto tempPlasticStrain = status->giveTempPlasticLatticeStrain();

        double le = static_cast < LatticeStructuralElement * > ( gp->giveElement() )->giveLength();

        double equivalentStrain = sqrt( pow(tempPlasticStrain.at(4), 2) + pow(tempPlasticStrain.at(5), 2.) + pow(tempPlasticStrain.at(6), 2) ) - wu / le;


        if ( equivalentStrain > status->giveKappaD() ) {
            tempKappaD = equivalentStrain;
        } else {
            tempKappaD = status->giveKappaD();
        }


        if ( tempKappaD  <= 0.  ) {
            omega = 0.;
        } else if ( tempKappaD > 0. ) {
            omega = 1. - exp( -tempKappaD * le / ( wf ) );
        } else {
            printf("Should not be here\n");
        }

        if ( omega > 1.0 ) {
            omega = 1.;
        } else if ( omega < 0.0 ) {
            omega = 0.;
        }

        stress *= ( 1. - omega );
    }

    status->letTempLatticeStrainBe(originalStrain);
    status->letTempReducedLatticeStrainBe(strain);
    status->letTempLatticeStressBe(stress);
    status->setTempKappaD(tempKappaD);
    status->setTempDamage(omega);

    return stress;
}


FloatMatrixF < 6, 6 >
LatticeFrameConcretePlastic::give3dFrameStiffnessMatrix(MatResponseMode rmode, GaussPoint * gp, TimeStep * atTime) const
{
    static_cast < LatticeFrameConcretePlasticStatus * > ( this->giveStatus(gp) );

    double g = this->e / ( 2. * ( 1. + this->nu ) );
    const double area       = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveArea(gp);
    const double shearareay = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveShearArea1(gp);
    const double shearareaz = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveShearArea2(gp);
    const double ik         = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveJ(gp);
    const double iy         = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveI1(gp);
    const double iz         = ( static_cast < LatticeStructuralElement * > ( gp->giveElement() ) )->giveI2(gp);


    FloatArrayF < 6 > d = {
        this->e * area,
        g * shearareay,
        g * shearareaz,
        g * ik,
        this->e * iy,
        this->e * iz
    };

    return diag(d);
}

LatticeFrameConcretePlasticStatus::LatticeFrameConcretePlasticStatus(int n, Domain *d, GaussPoint *g) :
    LatticeMaterialStatus(g)
{}
void
LatticeFrameConcretePlasticStatus::initTempStatus()
{
    LatticeMaterialStatus::initTempStatus();
    this->tempKappaD = this->kappaD;
    this->tempDamage = this->damage;
}

void
LatticeFrameConcretePlasticStatus::updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    LatticeMaterialStatus::updateYourself(atTime);
    this->kappaD = this->tempKappaD;
    this->damage = this->tempDamage;
}

void
LatticeFrameConcretePlasticStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    LatticeMaterialStatus::printOutputAt(file, tStep);

    fprintf(file, "plasticStrains ");
    for ( double s : this->plasticLatticeStrain ) {
        fprintf(file, "% .8e ", s);
    }
    fprintf(file, ", kappaD %.8e, damage %.8e \n", this->kappaD, this->damage);
    fprintf(file, "\n");
}

void
LatticeFrameConcretePlasticStatus::saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    LatticeMaterialStatus::saveContext(stream, mode);



    if ( !stream.write(& kappaD, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}



void
LatticeFrameConcretePlasticStatus::restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    LatticeMaterialStatus::restoreContext(stream, mode);

    if ( !stream.read(& kappaD, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}
}
