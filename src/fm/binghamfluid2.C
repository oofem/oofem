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

#include "binghamfluid2.h"
#include "fluiddynamicmaterial.h"
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
static bool __dummy_BinghamFluidMaterial2_alt = GiveClassFactory().registerMaterial("binghamfluid2", matCreator< BinghamFluidMaterial2 > );

BinghamFluidMaterial2 :: BinghamFluidMaterial2(int n, Domain *d) : FluidDynamicMaterial(n, d),
    mu_0(0.),
    tau_0(0.),
    tau_c(0.),
    mu_inf(1.e6),
    stressGrowthRate(BINGHAM_DEFAULT_STRESS_GROWTH_RATE)
{ }


IRResultType
BinghamFluidMaterial2 :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

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
    return FluidDynamicMaterial :: initializeFrom(ir);
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
BinghamFluidMaterial2 :: giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep)
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
BinghamFluidMaterial2 :: give(int aProperty, GaussPoint *gp)
{
    if ( aProperty == Viscosity ) {
        return mu_0;
    } else if ( aProperty == YieldStress ) {
        return tau_0;
    } else {
        return FluidDynamicMaterial :: give(aProperty, gp);
    }
}


MaterialStatus *
BinghamFluidMaterial2 :: CreateStatus(GaussPoint *gp) const
{
    return new BinghamFluidMaterial2Status(1, this->giveDomain(), gp);
}


void
BinghamFluidMaterial2 :: computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    int size = eps.giveSize();
    answer.resize(size);
    BinghamFluidMaterial2Status *status = static_cast< BinghamFluidMaterial2Status * >( this->giveStatus(gp) );
    MaterialMode mmode = gp->giveMaterialMode();
    FloatArray epsd;
    double gamma, tau, _nu;

    // determine actual viscosity
    this->computeDeviatoricStrain(epsd, eps, mmode);
    // determine shear strain magnitude
    gamma = this->computeDevStrainMagnitude(mmode, epsd);

#ifdef BINGHAM_ALT
    _nu = computeActualViscosity(tau_0, gamma);
    this->computeDeviatoricStress(answer, epsd, _nu, mmode);
    tau = this->computeDevStressMagnitude(mmode, answer);

    //printf ("_nu %e gamma %e\n", _nu, gamma);
#else
    // compute trial state
    this->computeDeviatoricStress(answer, epsd, this->mu_inf, mmode);
    // check if state allowed
    tau = this->computeDevStressMagnitude(mmode, answer);
    if ( tau > this->tau_c ) {
        _nu = this->computeActualViscosity(tau, gamma);
        this->computeDeviatoricStress(answer, epsd, _nu, mmode);
        tau = this->computeDevStressMagnitude(mmode, answer);
    }

#endif
    // update status
    status->letTempDeviatoricStrainVectorBe(epsd);
    status->letDeviatoricStressVectorBe(answer);
    status->letTempDevStrainMagnitudeBe(gamma);
    status->letTempDevStressMagnitudeBe(tau);
}

void
BinghamFluidMaterial2 :: giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp,
                                                       TimeStep *tStep)
{
    BinghamFluidMaterial2Status *status = static_cast< BinghamFluidMaterial2Status * >( this->giveStatus(gp) );
    MaterialMode mmode = gp->giveMaterialMode();
    const FloatArray &epsd = status->giveTempDeviatoricStrainVector(); //status->giveTempDeviatoricStrainVector();
    ///@note This variable was actually never used:
    //double tau = status->giveTempDevStressMagnitude();
    double gamma = status->giveTempDevStrainMagnitude(); //status->giveTempDevStrainMagnitude();
    // determine actual viscosity
    double gamma2 = gamma * gamma;

    if ( mmode == _2dFlow ) {
        answer.resize(3, 3);
        answer.zero();

        double dgde1, dgde2, dgde3;
        double dmudg, mu;

#if 0
        double _nu = computeActualViscosity(tau_0, gamma);

        answer.at(1, 1) = answer.at(2, 2) = 2.0 * _nu;
        answer.at(3, 3) = _nu;

        //answer.at(1,1) = answer.at(2,2) = answer.at(3,3) = 2.0*_nu;
        //answer.at(4,4) = _nu;
        return;
#endif

        if ( ( mode == ElasticStiffness ) || ( mode == SecantStiffness ) ) {
            double _nu = computeActualViscosity(tau_0, gamma);
            answer.at(1, 1) = answer.at(2, 2) = 2.0 * _nu;
            answer.at(3, 3) = _nu;
            return;
        } else { // tangent stiffness
            if ( gamma < BINGHAM_MIN_SHEAR_RATE ) {
                dmudg = dgde1 = dgde2 = dgde3 = 0.0;
                mu = computeActualViscosity(tau_0, gamma);
            } else {
                dmudg = ( -1.0 ) * tau_0 * ( 1.0 - exp(-this->stressGrowthRate * gamma) ) / gamma2 +
                        tau_0 *this->stressGrowthRate *exp(-this->stressGrowthRate *gamma) / gamma;
                mu = mu_0 + tau_0 * ( 1. - exp(-this->stressGrowthRate * gamma) ) / gamma;

#if 1
                //dgde1 = 2.0 * fabs( epsd.at(1) ) / gamma;
                //dgde2 = 2.0 * fabs( epsd.at(2) ) / gamma;
                //dgde3 = 1.0 * fabs( epsd.at(3) ) / gamma;
                dgde1 = 0.5 * fabs( epsd.at(1) ) / gamma;
                dgde2 = 0.5 * fabs( epsd.at(2) ) / gamma;
                dgde3 = 1.0 * fabs( epsd.at(3) ) / gamma;
#else
                //dgde1 = 2.0 * epsd.at(1) / gamma;
                //dgde2 = 2.0 * epsd.at(2) / gamma;
                //dgde3 = 1.0 * epsd.at(3) / gamma;
#endif
            }

            answer.at(1, 1) = 2.0 * epsd.at(1) * dmudg * dgde1 + 2.0 * mu;
            answer.at(1, 2) = 2.0 * epsd.at(1) * dmudg * dgde2;
            answer.at(1, 3) = 2.0 * epsd.at(1) * dmudg * dgde3;

            answer.at(2, 1) = 2.0 * epsd.at(2) * dmudg * dgde1;
            answer.at(2, 2) = 2.0 * epsd.at(2) * dmudg * dgde2 + 2.0 * mu;
            answer.at(2, 3) = 2.0 * epsd.at(2) * dmudg * dgde3;

            answer.at(3, 1) = epsd.at(3) * dmudg * dgde1;
            answer.at(3, 2) = epsd.at(3) * dmudg * dgde2;
            answer.at(3, 3) = epsd.at(3) * dmudg * dgde3 + mu;

            return;
        }
    } else if ( mmode == _2dAxiFlow ) {
        answer.resize(4, 4);
        answer.zero();


        double dgde1, dgde2, dgde3, dgde4;
        double dmudg, mu;

        if ( 0 ) {
            double _nu = computeActualViscosity(tau_0, gamma);

            answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 2.0 * _nu;
            answer.at(4, 4) = _nu;

            //answer.at(1,1) = answer.at(2,2) = answer.at(3,3) = 2.0*_nu;
            //answer.at(4,4) = _nu;
            return;
        }

        if ( ( mode == ElasticStiffness ) || ( mode == SecantStiffness ) ) {
            double _nu = computeActualViscosity(tau_0, gamma);
            answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 2.0 * _nu;
            answer.at(4, 4) = _nu;
            return;
        } else { // tangent stiffness
            if ( gamma < BINGHAM_MIN_SHEAR_RATE ) {
                dmudg = dgde1 = dgde2 = dgde3 = dgde4 = 0.0;
                mu = computeActualViscosity(tau_0, gamma);
            } else {
                dmudg = ( -1.0 ) * tau_0 * ( 1.0 - exp(-this->stressGrowthRate * gamma) ) / gamma2 +
                        tau_0 *this->stressGrowthRate *exp(-this->stressGrowthRate *gamma) / gamma;
                mu = mu_0 + tau_0 * ( 1. - exp(-this->stressGrowthRate * gamma) ) / gamma;

#if 1
                dgde1 = 2.0 * fabs( epsd.at(1) ) / gamma;
                dgde2 = 2.0 * fabs( epsd.at(2) ) / gamma;
                dgde3 = 2.0 * fabs( epsd.at(3) ) / gamma;
                dgde4 = 1.0 * fabs( epsd.at(4) ) / gamma;
#else
                dgde1 = 2.0 * epsd.at(1) / gamma;
                dgde2 = 2.0 * epsd.at(2) / gamma;
                dgde3 = 2.0 * epsd.at(3) / gamma;
                dgde4 = 1.0 * epsd.at(4) / gamma;
#endif
            }

            answer.at(1, 1) = 2.0 * epsd.at(1) * dmudg * dgde1 + 2.0 * mu;
            answer.at(1, 2) = 2.0 * epsd.at(1) * dmudg * dgde2;
            answer.at(1, 3) = 2.0 * epsd.at(1) * dmudg * dgde3;
            answer.at(1, 4) = 2.0 * epsd.at(1) * dmudg * dgde4;

            answer.at(2, 1) = 2.0 * epsd.at(2) * dmudg * dgde1;
            answer.at(2, 2) = 2.0 * epsd.at(2) * dmudg * dgde2 + 2.0 * mu;
            answer.at(2, 3) = 2.0 * epsd.at(2) * dmudg * dgde3;
            answer.at(2, 4) = 2.0 * epsd.at(2) * dmudg * dgde4;

            answer.at(3, 1) = 2.0 * epsd.at(3) * dmudg * dgde1;
            answer.at(3, 2) = 2.0 * epsd.at(3) * dmudg * dgde2;
            answer.at(3, 3) = 2.0 * epsd.at(3) * dmudg * dgde3 + 2.0 * mu;
            answer.at(3, 4) = 2.0 * epsd.at(3) * dmudg * dgde4;


            answer.at(4, 1) = epsd.at(4) * dmudg * dgde1;
            answer.at(4, 2) = epsd.at(4) * dmudg * dgde2;
            answer.at(4, 3) = epsd.at(4) * dmudg * dgde3;
            answer.at(4, 4) = epsd.at(4) * dmudg * dgde4 + mu;

            return;
        }
    } else if ( gp->giveMaterialMode() == _3dFlow ) {
        answer.resize(6, 6);
        answer.zero();

        FloatArray dgde(6);
        double dmudg, mu;

        if ( gamma < BINGHAM_MIN_SHEAR_RATE ) {
            dmudg = 0.0;
            dgde.zero();
            mu = computeActualViscosity(tau_0, gamma);
        } else {
            dmudg = ( -1.0 ) * tau_0 * ( 1.0 - exp(-this->stressGrowthRate * gamma) ) / gamma2 +
                    tau_0 *this->stressGrowthRate *exp(-this->stressGrowthRate *gamma) / gamma;
            mu = mu_0 + tau_0 * ( 1. - exp(-this->stressGrowthRate * gamma) ) / gamma;

            dgde.at(1) = 2.0 * epsd.at(1) / gamma;
            dgde.at(2) = 2.0 * epsd.at(2) / gamma;
            dgde.at(3) = 2.0 * epsd.at(3) / gamma;
            dgde.at(4) = 1.0 * epsd.at(4) / gamma;
            dgde.at(5) = 1.0 * epsd.at(5) / gamma;
            dgde.at(6) = 1.0 * epsd.at(6) / gamma;
        }

        for ( int i = 1; i <= 6; i++ ) {
            answer.at(1, i) = fabs( 2.0 * epsd.at(1) * dmudg * dgde.at(i) );
            answer.at(2, i) = fabs( 2.0 * epsd.at(2) * dmudg * dgde.at(i) );
            answer.at(3, i) = fabs( 2.0 * epsd.at(3) * dmudg * dgde.at(i) );
            answer.at(4, i) = fabs( epsd.at(4) * dmudg * dgde.at(i) );
            answer.at(5, i) = fabs( epsd.at(5) * dmudg * dgde.at(i) );
            answer.at(6, i) = fabs( epsd.at(6) * dmudg * dgde.at(i) );
        }

        answer.at(1, 1) += 2.0 * mu;
        answer.at(2, 2) += 2.0 * mu;
        answer.at(3, 3) += 2.0 * mu;
        answer.at(4, 4) += mu;
        answer.at(5, 5) += mu;
        answer.at(6, 6) += mu;

        return;
    }  else {
        OOFEM_ERROR("unsupportted material mode");
    }
}


int
BinghamFluidMaterial2 :: checkConsistency()
{
    if ( domain->giveEngngModel()->giveEquationScalingFlag() ) {
        double scale;
        scale = domain->giveEngngModel()->giveVariableScale(VST_Density);
        propertyDictionary.at('d') /= scale;

        scale = domain->giveEngngModel()->giveVariableScale(VST_Viscosity);
        this->mu_0 /= scale;
        this->tau_0 /= scale;
    }

    return 1;
}

double
BinghamFluidMaterial2 :: computeActualViscosity(double Tau, double shearRate)
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
    if ( Tau <= tau_c ) {
        return this->mu_inf;
    } else {
        return mu_0 + ( tau_0 / shearRate );
    }

#endif
}


double
BinghamFluidMaterial2 :: computeDevStrainMagnitude(MaterialMode mmode, const FloatArray &epsd)
{
    double _val = 0.0;
    if ( mmode == _2dFlow ) {
        // _val = 2.0 * ( epsd.at(1) * epsd.at(1) + epsd.at(2) * epsd.at(2) ) + epsd.at(3) * epsd.at(3);
        _val = 0.5 * ( epsd.at(1) * epsd.at(1) + epsd.at(2) * epsd.at(2) ) + epsd.at(3) * epsd.at(3);
    } else if ( mmode == _2dAxiFlow ) {
        _val = 2.0 * ( epsd.at(1) * epsd.at(1) + epsd.at(2) * epsd.at(2) +
                      epsd.at(3) * epsd.at(3) ) + epsd.at(4) * epsd.at(4);
    } else if ( mmode == _3dFlow ) {
        _val = 2.0 * ( epsd.at(1) * epsd.at(1) + epsd.at(2) * epsd.at(2) + epsd.at(3) * epsd.at(3) )
               + epsd.at(4) * epsd.at(4) + epsd.at(5) * epsd.at(5) + epsd.at(6) * epsd.at(6);
    } else {
        OOFEM_ERROR("unsupported material mode");
    }

    return sqrt(_val);
}

double
BinghamFluidMaterial2 :: computeDevStressMagnitude(MaterialMode mmode, const FloatArray &sigd)
{
    double _val = 0.0;
    if ( mmode == _2dFlow ) {
        _val = 0.5 * ( sigd.at(1) * sigd.at(1) + sigd.at(2) * sigd.at(2) + 2.0 * sigd.at(3) * sigd.at(3) );
    } else if ( mmode == _2dAxiFlow ) {
        _val = 0.5 * ( sigd.at(1) * sigd.at(1) +
                      sigd.at(2) * sigd.at(2) +
                      sigd.at(3) * sigd.at(3) +
                      2.0 * sigd.at(4) * sigd.at(4) );
    } else if ( mmode == _3dFlow ) {
        _val = 0.5 * ( sigd.at(1) * sigd.at(1) + sigd.at(2) * sigd.at(2) + sigd.at(3) * sigd.at(3) +
                      2.0 * sigd.at(4) * sigd.at(4) + 2.0 * sigd.at(5) * sigd.at(5) + 2.0 * sigd.at(6) * sigd.at(6) );
    } else {
        OOFEM_ERROR("unsupported material mode");
    }

    return sqrt(_val);
}

void
BinghamFluidMaterial2 :: computeDeviatoricStrain(FloatArray &answer, const FloatArray &eps, MaterialMode mmode)
{
    if ( mmode ==  _2dFlow ) {
        //double ekk=(eps.at(1)+eps.at(2))/3.0;
        double ekk = 0.0;

        answer = {
            eps.at(1) - ekk,
            eps.at(2) - ekk,
            eps.at(3)
        };
    } else if ( mmode == _2dAxiFlow ) {
        //double ekk=(eps.at(1)+eps.at(2)+eps.at(3))/3.0;
        double ekk = 0.0;

        answer = {
            eps.at(1) - ekk,
            eps.at(2) - ekk,
            eps.at(3) - ekk,
            eps.at(4)
        };
    } else if ( mmode == _3dFlow ) {
        //double ekk=(eps.at(1)+eps.at(2)+eps.at(3))/3.0;
        double ekk = 0.0;

        answer  = {
            eps.at(1) - ekk,
            eps.at(2) - ekk,
            eps.at(3) - ekk,
            eps.at(4),
            eps.at(5),
            eps.at(6)
        };
    } else {
        OOFEM_ERROR("unsupported material mode");
    }
}


void
BinghamFluidMaterial2 :: computeDeviatoricStress(FloatArray &answer, const FloatArray &deps,
                                                 double _nu, MaterialMode mmode)
{
    if ( mmode == _2dFlow ) {
        answer = {
            2.0 * _nu * ( deps.at(1) ),
            2.0 * _nu * ( deps.at(2) ),
            deps.at(3) * _nu
        };
    } else if ( mmode == _2dAxiFlow ) {
        answer = {
            answer.at(1) = 2.0 * _nu * ( deps.at(1) ),
            answer.at(2) = 2.0 * _nu * ( deps.at(2) ),
            answer.at(3) = 2.0 * _nu * ( deps.at(3) ),
            answer.at(4) = deps.at(4) * _nu,
        };
    } else if ( mmode == _3dFlow ) {
        answer = {
            2.0 * _nu * ( deps.at(1) ),
            2.0 * _nu * ( deps.at(2) ),
            2.0 * _nu * ( deps.at(3) ),
            deps.at(4) * _nu,
            deps.at(5) * _nu,
            deps.at(6) * _nu
        };
    } else {
        OOFEM_ERROR("unsupported material mode");
    }
}


BinghamFluidMaterial2Status :: BinghamFluidMaterial2Status(int n, Domain *d, GaussPoint *g) :
    FluidDynamicMaterialStatus(n, d, g)
{
    MaterialMode mmode = gp->giveMaterialMode();
    int _size = 0;

    devStrainMagnitude = temp_devStrainMagnitude = 0.0;
    devStressMagnitude = temp_devStressMagnitude = 0.0;

    if ( mmode == _2dFlow ) {
        _size = 3;
    } else if ( mmode == _2dAxiFlow ) {
        _size = 4;
    } else if ( mmode == _3dFlow ) {
        _size = 6;
    } else {
        OOFEM_ERROR("unsupported material mode");
    }

    deviatoricStrainRateVector.resize(_size);
    deviatoricStrainRateVector.zero();
    deviatoricStressVector.resize(_size);
    deviatoricStressVector.zero();
    temp_deviatoricStrainVector = deviatoricStrainRateVector;
}

void
BinghamFluidMaterial2Status :: printOutputAt(FILE *File, TimeStep *tStep)
// Prints the strains and stresses on the data file.
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
// Performs end-of-step updates.
{
    FluidDynamicMaterialStatus :: updateYourself(tStep);

    devStrainMagnitude = temp_devStrainMagnitude;
    devStressMagnitude = temp_devStressMagnitude;
    deviatoricStrainRateVector = temp_deviatoricStrainVector;
}


void
BinghamFluidMaterial2Status :: initTempStatus()
//
// initialize record at the begining of new load step
//
{
    FluidDynamicMaterialStatus :: initTempStatus();

    temp_devStrainMagnitude = devStrainMagnitude;
    temp_devStressMagnitude = devStressMagnitude;
    temp_deviatoricStrainVector = deviatoricStrainRateVector;
}


contextIOResultType
BinghamFluidMaterial2Status :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full ms context (saves state variables, that completely describe
// current state)
// saving the data in dictionary is left to material (yield crit. level).
{
    contextIOResultType iores;

    if ( ( iores = FluidDynamicMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(devStrainMagnitude) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(devStressMagnitude) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


contextIOResultType
BinghamFluidMaterial2Status :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = FluidDynamicMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(devStrainMagnitude) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(devStressMagnitude) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
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
