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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "fm/Materials/binghamfluid2.h"
#include "fm/Materials/fluiddynamicmaterial.h"
#include "domain.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "engngm.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"

#include <cstdlib>

namespace oofem {
#define BINGHAM_ALT 1
#define BINGHAM_MIN_SHEAR_RATE     1.e-10

REGISTER_Material(BinghamFluidMaterial2);
///@todo Remove the alternative ID. Just stick to "binghamfluid".
// static bool __dummy_BinghamFluidMaterial2_alt __attribute__((unused)) = GiveClassFactory().registerMaterial("binghamfluid2", matCreator< BinghamFluidMaterial2 > );
REGISTER_Material_Alt(BinghamFluidMaterial2, binghamfluid2);

BinghamFluidMaterial2 :: BinghamFluidMaterial2(int n, Domain *d) : FluidDynamicMaterial(n, d)
{ }


void
BinghamFluidMaterial2 :: initializeFrom(InputRecord &ir)
{
    FluidDynamicMaterial :: initializeFrom(ir);
    // we use rather object's member data than to store data into slow
    // key-val dictionary with lot of memory allocations
    IR_GIVE_FIELD(ir, mu_0, _IFT_BinghamFluidMaterial2_mu0);
    IR_GIVE_FIELD(ir, tau_0, _IFT_BinghamFluidMaterial2_tau0);
    mu_inf = 1.e6;
    IR_GIVE_OPTIONAL_FIELD(ir, mu_inf, _IFT_BinghamFluidMaterial2_muinf);
    stressGrowthRate = BINGHAM_DEFAULT_STRESS_GROWTH_RATE;
    IR_GIVE_OPTIONAL_FIELD(ir, stressGrowthRate, _IFT_BinghamFluidMaterial2_stressGrowthRate);
    tau_c = tau_0 * mu_inf / ( mu_inf - mu_0 );
    //tau_c = tau_0;
}


void
BinghamFluidMaterial2 :: giveInputRecord(DynamicInputRecord &input)
{
    FluidDynamicMaterial :: giveInputRecord(input);
    input.setField(this->mu_0, _IFT_BinghamFluidMaterial2_mu0);
    input.setField(this->tau_0, _IFT_BinghamFluidMaterial2_tau0);
    input.setField(this->mu_inf, _IFT_BinghamFluidMaterial2_muinf);
    input.setField(this->stressGrowthRate, _IFT_BinghamFluidMaterial2_stressGrowthRate);
}

double
BinghamFluidMaterial2 :: giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep) const
{
    BinghamFluidMaterial2Status *status = static_cast< BinghamFluidMaterial2Status * >( this->giveStatus(gp) );
    //double temp_tau=status->giveTempDevStressMagnitude();
    //double temp_gamma=status->giveTempDevStrainMagnitude();
    //determine actual viscosity
    //return this->computeActualViscosity(temp_tau, temp_gamma);
#ifdef BINGHAM_ALT
    double gamma = status->giveTempDevStrainMagnitude(); //status->giveTempDevStrainMagnitude();
    return computeActualViscosity(tau_0, gamma);

 #if 0
    const FloatArray &epsd = status->giveTempDeviatoricStrainVector(); //status->giveTempDeviatoricStrainVector();
    double gamma2 = gamma * gamma;
    double dmudg, dgde1, dgde2, dgde3, mu;
    if ( gamma < BINGHAM_MIN_SHEAR_RATE ) {
        dmudg = dgde1 = dgde2 = dgde3 = 0.0;
        mu = computeActualViscosity(tau_0, gamma);
    } else {
        dmudg = ( -1.0 ) * tau_0 * ( 1.0 - exp(-this->stressGrowthRate * gamma) ) / gamma2 +
                tau_0 *this->stressGrowthRate *exp(-this->stressGrowthRate *gamma) / gamma;
        mu = mu_0 + tau_0 * ( 1. - exp(-this->stressGrowthRate * gamma) ) / gamma;

        dgde1 = 2.0 * fabs( epsd.at(1) ) / gamma;
        dgde2 = 2.0 * fabs( epsd.at(2) ) / gamma;
        dgde3 = 1.0 * fabs( epsd.at(3) ) / gamma;
    }

    return min( min( ( epsd.at(1) * dmudg * dgde1 + mu ), ( epsd.at(2) * dmudg * dgde2 + mu ) ),
               ( epsd.at(3) * dmudg * dgde3 + mu ) );

 #endif

#else
    if ( temp_tau < tau_c ) {
        return mu_inf;
    } else {
        return mu_0;
    }

#endif
}


double
BinghamFluidMaterial2 :: give(int aProperty, GaussPoint *gp) const
{
    if ( aProperty == Viscosity ) {
        return mu_0;
    } else if ( aProperty == YieldStress ) {
        return tau_0;
    } else {
        return FluidDynamicMaterial :: give(aProperty, gp);
    }
}


std::unique_ptr<MaterialStatus> 
BinghamFluidMaterial2 :: CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<BinghamFluidMaterial2Status>(gp);
}


FloatArrayF<6>
BinghamFluidMaterial2 :: computeDeviatoricStress3D(const FloatArrayF<6> &eps, GaussPoint *gp, TimeStep *tStep) const
{
    BinghamFluidMaterial2Status *status = static_cast< BinghamFluidMaterial2Status * >( this->giveStatus(gp) );

    // determine actual viscosity
    auto epsd = this->computeDeviatoricStrain(eps);
    // determine shear strain magnitude
    double gamma = this->computeDevStrainMagnitude(epsd);

#ifdef BINGHAM_ALT
    double nu = computeActualViscosity(tau_0, gamma);
    auto stress = this->computeDeviatoricStress(epsd, nu);
    double tau = this->computeDevStressMagnitude(stress);

    //printf ("nu %e gamma %e\n", nu, gamma);
#else
    // compute trial state
    auto stress = this->computeDeviatoricStress(epsd, this->mu_inf);
    // check if state allowed
    double tau = this->computeDevStressMagnitude(stress);
    if ( tau > this->tau_c ) {
        double nu = this->computeActualViscosity(tau, gamma);
        auto stress = this->computeDeviatoricStress(stress, epsd, nu);
        tau = this->computeDevStressMagnitude(stress);
    }

#endif
    // update status
    status->letTempDeviatoricStrainVectorBe(epsd);
    status->letDeviatoricStressVectorBe(stress);
    status->letTempDevStrainMagnitudeBe(gamma);
    status->letTempDevStressMagnitudeBe(tau);

    return stress;
}

FloatMatrixF<6,6>
BinghamFluidMaterial2 :: computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    BinghamFluidMaterial2Status *status = static_cast< BinghamFluidMaterial2Status * >( this->giveStatus(gp) );
    const auto &epsd = status->giveTempDeviatoricStrainVector();
    double gamma = status->giveTempDevStrainMagnitude();

    FloatMatrixF<6,6> d;
    double mu;
    if ( gamma < BINGHAM_MIN_SHEAR_RATE ) {
        mu = computeActualViscosity(tau_0, gamma);
    } else {
        double dmudg = - tau_0 * ( 1.0 - exp(-this->stressGrowthRate * gamma) ) / (gamma * gamma) +
                       tau_0 * this->stressGrowthRate * exp(-this->stressGrowthRate * gamma) / gamma;
        mu = mu_0 + tau_0 * ( 1. - exp(-this->stressGrowthRate * gamma) ) / gamma;

        FloatArrayF<6> dgde = {
            2.0 * epsd.at(1) / gamma,
            2.0 * epsd.at(2) / gamma,
            2.0 * epsd.at(3) / gamma,
            1.0 * epsd.at(4) / gamma,
            1.0 * epsd.at(5) / gamma,
            1.0 * epsd.at(6) / gamma,
        };

        for ( int i = 1; i <= 6; i++ ) {
            d.at(1, i) = std::abs( 2.0 * epsd.at(1) * dmudg * dgde.at(i) );
            d.at(2, i) = std::abs( 2.0 * epsd.at(2) * dmudg * dgde.at(i) );
            d.at(3, i) = std::abs( 2.0 * epsd.at(3) * dmudg * dgde.at(i) );
            d.at(4, i) = std::abs( epsd.at(4) * dmudg * dgde.at(i) );
            d.at(5, i) = std::abs( epsd.at(5) * dmudg * dgde.at(i) );
            d.at(6, i) = std::abs( epsd.at(6) * dmudg * dgde.at(i) );
        }
    }

    d.at(1, 1) += 2.0 * mu;
    d.at(2, 2) += 2.0 * mu;
    d.at(3, 3) += 2.0 * mu;
    d.at(4, 4) += mu;
    d.at(5, 5) += mu;
    d.at(6, 6) += mu;
    return d;
}


int
BinghamFluidMaterial2 :: checkConsistency()
{
    if ( domain->giveEngngModel()->giveEquationScalingFlag() ) {
        double scale = domain->giveEngngModel()->giveVariableScale(VST_Density);
        propertyDictionary.at('d') /= scale;

        scale = domain->giveEngngModel()->giveVariableScale(VST_Viscosity);
        this->mu_0 /= scale;
        this->tau_0 /= scale;
    }

    return 1;
}

double
BinghamFluidMaterial2 :: computeActualViscosity(double tau, double shearRate) const
{
#ifdef BINGHAM_ALT
    if ( tau_0 > 0.0 ) {
        shearRate = max(shearRate, BINGHAM_MIN_SHEAR_RATE);
        return ( mu_0 + tau_0 * ( 1. - exp(-this->stressGrowthRate * shearRate) ) / shearRate );
    } else {
        // newtonian flow
        return mu_0;
    }

#else
    if ( tau <= tau_c ) {
        return this->mu_inf;
    } else {
        return mu_0 + ( tau_0 / shearRate );
    }

#endif
}


double
BinghamFluidMaterial2 :: computeDevStrainMagnitude(const FloatArrayF<6> &epsd)
{
    double val = 2.0 * ( epsd[0] * epsd[0] + epsd[1] * epsd[1] + epsd[2] * epsd[2] ) + 
                epsd[3] * epsd[3] + epsd[4] * epsd[4] + epsd[5] * epsd[5];
    return sqrt(val);
}

double
BinghamFluidMaterial2 :: computeDevStressMagnitude(const FloatArrayF<6> &sigd)
{
    double val = 0.5 * ( sigd[0] * sigd[0] + sigd[1] * sigd[1] + sigd[2] * sigd[2] +
                2.0 * (sigd[3] * sigd[3] + sigd[4] * sigd[4] + sigd[5] * sigd[5]) );
    return sqrt(val);
}

FloatArrayF<6>
BinghamFluidMaterial2 :: computeDeviatoricStrain(const FloatArrayF<6> &eps)
{
    //double ekk=(eps.at(1)+eps.at(2)+eps.at(3))/3.0;
    double ekk = 0.0;
    return {
        eps[0] - ekk,
        eps[1] - ekk,
        eps[2] - ekk,
        eps[3],
        eps[4],
        eps[5]
    };
}


FloatArrayF<6>
BinghamFluidMaterial2 :: computeDeviatoricStress(const FloatArrayF<6> &deps, double nu)
{
    return {
        2.0 * nu * ( deps[0] ),
        2.0 * nu * ( deps[1] ),
        2.0 * nu * ( deps[2] ),
        deps[3] * nu,
        deps[4] * nu,
        deps[5] * nu
    };
}


BinghamFluidMaterial2Status :: BinghamFluidMaterial2Status(GaussPoint *g) :
    FluidDynamicMaterialStatus(g)
{}

void
BinghamFluidMaterial2Status :: printOutputAt(FILE *File, TimeStep *tStep) const
{
    fprintf(File, " strains ");
    for ( double e: deviatoricStrainRateVector ) {
        fprintf( File, " %.4e", e );
    }

    fprintf(File, "\n deviatoric stresses");
    for ( double e: deviatoricStressVector ) {
        fprintf( File, " %.4e", e );
    }

    fprintf(File, "\n          status { gamma %e, tau %e }", devStrainMagnitude, devStressMagnitude);

    fprintf(File, "\n");
}

void
BinghamFluidMaterial2Status :: updateYourself(TimeStep *tStep)
{
    FluidDynamicMaterialStatus :: updateYourself(tStep);

    devStrainMagnitude = temp_devStrainMagnitude;
    devStressMagnitude = temp_devStressMagnitude;
    deviatoricStrainRateVector = temp_deviatoricStrainVector;
}


void
BinghamFluidMaterial2Status :: initTempStatus()
{
    FluidDynamicMaterialStatus :: initTempStatus();

    temp_devStrainMagnitude = devStrainMagnitude;
    temp_devStressMagnitude = devStressMagnitude;
    temp_deviatoricStrainVector = deviatoricStrainRateVector;
}


void
BinghamFluidMaterial2Status :: saveContext(DataStream &stream, ContextMode mode)
{
    FluidDynamicMaterialStatus :: saveContext(stream, mode);

    if ( !stream.write(devStrainMagnitude) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(devStressMagnitude) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


void
BinghamFluidMaterial2Status :: restoreContext(DataStream &stream, ContextMode mode)
{
    FluidDynamicMaterialStatus :: restoreContext(stream, mode);

    if ( !stream.read(devStrainMagnitude) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(devStressMagnitude) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


#if 0
void
BinghamFluidMaterial2 :: __debug(GaussPoint *gp, TimeStep *tStep)
{
    BinghamFluidMaterial2Status *status = static_cast< BinghamFluidMaterial2Status * >( this->giveStatus(gp) );
    const FloatArray &epsd = status->giveTempDeviatoricStrainVector();
    const FloatArray &sigd = status->giveTempDeviatoricStrainVector()
    for ( int i = 1; i <= nincr; i++ ) {
        eps.add(eps_i);
        computeDeviatoricStressVector(tau, gp, eps, tStep);
        giveDeviatoricStiffnessMatrix(d, TangentStiffness, gp, tStep);
        tau_t.beProductOf(d, eps_i);
        tau_t.add(tau_p);
        //tau.printYourself();
        //tau_t.printYourself();
        //d.printYourself();
        printf( "%e %e %e  %e %e %e %e %e %e\n", eps.at(1), eps.at(2), eps.at(3), tau.at(1), tau.at(2), tau.at(3), tau_t.at(1), tau_t.at(2), tau_t.at(3) );
        tau_p = tau_t;
    }
}
#endif
} // end namespace oofem
