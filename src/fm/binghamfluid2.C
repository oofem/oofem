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

#include "binghamfluid2.h"
#include "fluiddynamicmaterial.h"
#include "domain.h"
#include "flotmtrx.h"
#include "gausspnt.h"
#include "engngm.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"


namespace oofem {
#define BINGHAM_ALT 1
#define BINGHAM_MIN_SHEAR_RATE     1.e-12

int
BinghamFluidMaterial2 :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    if ( ( mode == _2dFlow ) || ( mode == _3dFlow ) ) {
        return 1;
    }


    return 0;
}


IRResultType
BinghamFluidMaterial2 :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->FluidDynamicMaterial :: initializeFrom(ir);
    // we use rather object's member data than to store data into slow
    // key-val dictionary with lot of memory allocations

    IR_GIVE_FIELD(ir, mu_0, IFT_BinghamFluidMaterial_mu0, "mu0"); // Macro
    IR_GIVE_FIELD(ir, tau_0, IFT_BinghamFluidMaterial_tau0, "tau0"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, mu_inf, IFT_BinghamFluidMaterial_muinf, "muinf"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, stressGrowthRate, IFT_BinghamFluidMaterial_stressGrowthRate, "stressgrowthrate"); // Macro
    tau_c = tau_0 * mu_inf / ( mu_inf - mu_0 );
    //tau_c = tau_0;
    return IRRT_OK;
}


int
BinghamFluidMaterial2 :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    FluidDynamicMaterial :: giveInputRecordString(str, keyword);
    sprintf(buff, " mu0 %e tau0 %e", this->mu_0, this->tau_0);
    str += buff;

    return 1;
}


double
BinghamFluidMaterial2 :: giveCharacteristicValue(MatResponseMode mode,
        GaussPoint *gp,
        TimeStep *atTime)
{
    if ( mode == MRM_Density ) {
        return this->give('d', gp);
    } else if ( mode == MRM_Viscosity ) {
        BinghamFluidMaterial2Status *status = ( ( BinghamFluidMaterial2Status * ) this->giveStatus(gp) );
        //double temp_tau=status->giveTempDevStressMagnitude();
        //double temp_gamma=status->giveTempDevStrainMagnitude();
        //determine actual viscosity
        //return this->computeActualViscosity(temp_tau, temp_gamma);
#ifdef BINGHAM_ALT
        double gamma = status->giveTempDevStrainMagnitude(); //status->giveTempDevStrainMagnitude();
        return computeActualViscosity(tau_0, gamma);

        /*
         * const FloatArray &epsd = status->giveTempDeviatoricStrainVector(); //status->giveTempDeviatoricStrainVector();
         * double gamma2 = gamma * gamma;
         * double dmudg, dgde1, dgde2, dgde3, mu;
         * if ( gamma < BINGHAM_MIN_SHEAR_RATE ) {
         *  dmudg = dgde1 = dgde2 = dgde3 = 0.0;
         *  mu = computeActualViscosity(tau_0, gamma);
         * } else {
         *  dmudg = ( -1.0 ) * tau_0 * ( 1.0 - exp(-this->stressGrowthRate * gamma) ) / gamma2 +
         *          tau_0 *this->stressGrowthRate *exp(-this->stressGrowthRate * gamma) / gamma;
         *  mu = mu_0 + tau_0 * ( 1. - exp(-this->stressGrowthRate * gamma) ) / gamma;
         *
         *  dgde1 = 2.0 * fabs( epsd.at(1) ) / gamma;
         *  dgde2 = 2.0 * fabs( epsd.at(2) ) / gamma;
         *  dgde3 = 1.0 * fabs( epsd.at(3) ) / gamma;
         * }
         *
         * return min( min( ( epsd.at(1) * dmudg * dgde1 + mu ), ( epsd.at(2) * dmudg * dgde2 + mu ) ),
         *         ( epsd.at(3) * dmudg * dgde3 + mu ) );
         */

#else
        if ( temp_tau < tau_c ) {
            return mu_inf;
        } else {
            return mu_0;
        }

#endif
        //return this->mu_0;
    } else {
        return FluidDynamicMaterial :: giveCharacteristicValue(mode, gp, atTime);
    }
}


double
BinghamFluidMaterial2 :: give(int aProperty, GaussPoint *gp)
//
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
//
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
/*
 * creates new  material status  corresponding to this class
 */
{
    return new BinghamFluidMaterial2Status(1, this->giveDomain(), gp);
}


void
BinghamFluidMaterial2 :: computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    int size = eps.giveSize();
    answer.resize(size);
    BinghamFluidMaterial2Status *status = ( ( BinghamFluidMaterial2Status * ) this->giveStatus(gp) );
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
    status->letTempDeviatoricStressVectorBe(answer);
    status->letTempDevStrainMagnitudeBe(gamma);
    status->letTempDevStressMagnitudeBe(tau);
}

void
BinghamFluidMaterial2 :: giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp,
        TimeStep *atTime)
{
    BinghamFluidMaterial2Status *status = ( ( BinghamFluidMaterial2Status * ) this->giveStatus(gp) );
    MaterialMode mmode = gp->giveMaterialMode();
    const FloatArray &epsd = status->giveTempDeviatoricStrainVector();
    double gamma = status->giveTempDevStrainMagnitude();
    // determine actual viscosity
    double _nu;
    double gamma2 = gamma * gamma;
    double dmudg, mu;

    if ( mmode == _2dFlow ) {
        answer.resize(3, 3);
        answer.zero();

        double dgde1, dgde2, dgde3, dgde4;
        double dmudg, mu;

        if ( 0 ) {
            _nu = computeActualViscosity(tau_0, gamma);

            answer.at(1, 1) = answer.at(2, 2) = 2.0 * _nu;
            answer.at(3, 3) = _nu;

            //answer.at(1,1) = answer.at(2,2) = answer.at(3,3) = 2.0*_nu;
            //answer.at(4,4) = _nu;
            return;
        }

        if ( ( mode == ElasticStiffness ) || ( mode == SecantStiffness ) ) {
            _nu = computeActualViscosity(tau_0, gamma);
            answer.at(1, 1) = answer.at(2, 2) = 2.0 * _nu;
            answer.at(3, 3) = _nu;
            return;
        } else { // tangent stiffness

#ifdef BINGHAM_MIN_SHEAR_RATE
            if ( gamma < BINGHAM_MIN_SHEAR_RATE ) {
                dmudg = dgde1 = dgde2 = dgde3 = dgde4 = 0.0;
                mu = computeActualViscosity(tau_0, gamma);
            } else {
#endif
                dmudg = ( -1.0 ) * tau_0 * ( 1.0 - exp(-this->stressGrowthRate * gamma) ) / gamma2 +
                        tau_0 *this->stressGrowthRate *exp(-this->stressGrowthRate *gamma) / gamma;
                mu = mu_0 + tau_0 * ( 1. - exp(-this->stressGrowthRate * gamma) ) / gamma;

#if 1
                dgde1 = 2.0 * fabs( epsd.at(1) ) / gamma;
                dgde2 = 2.0 * fabs( epsd.at(2) ) / gamma;
                dgde3 = 1.0 * fabs( epsd.at(3) ) / gamma;
#else
                dgde1 = 2.0 * epsd.at(1) / gamma;
                dgde2 = 2.0 * epsd.at(2) / gamma;
                dgde3 = 1.0 * epsd.at(3) / gamma;
#endif
#ifdef BINGHAM_MIN_SHEAR_RATE
            }
#endif

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

        if ( 0 ) {
            _nu = computeActualViscosity(tau_0, gamma);

            answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 2.0 * _nu;
            answer.at(4, 4) = _nu;

            //answer.at(1,1) = answer.at(2,2) = answer.at(3,3) = 2.0*_nu;
            //answer.at(4,4) = _nu;
            return;
        }

        if ( ( mode == ElasticStiffness ) || ( mode == SecantStiffness ) ) {
            _nu = computeActualViscosity(tau_0, gamma);
            answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 2.0 * _nu;
            answer.at(4, 4) = _nu;
            return;
        } else { // tangent stiffness

#ifdef BINGHAM_MIN_SHEAR_RATE
            if ( gamma < BINGHAM_MIN_SHEAR_RATE ) {
                dmudg = dgde1 = dgde2 = dgde3 = dgde4 = 0.0;
                mu = computeActualViscosity(tau_0, BINGHAM_MIN_SHEAR_RATE);
            } else {
#endif
                dmudg = ( -1.0 ) * tau_0 * ( 1.0 - exp(-this->stressGrowthRate * gamma) ) / gamma2 +
                        tau_0 *this->stressGrowthRate *exp(-this->stressGrowthRate *gamma) / gamma;
                mu = mu_0 + tau_0 * ( 1. - exp(-this->stressGrowthRate * gamma) ) / gamma;

#if 1
                dgde1 = (1./6.)*2.0*fabs( epsd.at(1) );
                dgde2 = (1./6.)*2.0*fabs( epsd.at(2) );
                dgde3 = (1./6.)*2.0*fabs( epsd.at(3) );
                dgde4 = (1./6.)*2.0*fabs( epsd.at(4) );
#else
                dgde1 = 2.0 * epsd.at(1) / gamma;
                dgde2 = 2.0 * epsd.at(2) / gamma;
                dgde3 = 2.0 * epsd.at(3) / gamma;
                dgde4 = 1.0 * epsd.at(4) / gamma;
#endif
#ifdef BINGHAM_MIN_SHEAR_RATE
            }
#endif

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

        FloatArray dgde (6);

#ifdef BINGHAM_MIN_SHEAR_RATE
        if ( gamma < BINGHAM_MIN_SHEAR_RATE ) {
            dmudg = 0.0;
            dgde.zero();
            mu = computeActualViscosity(tau_0, gamma);
        } else {
#endif
            dmudg = ( -1.0 ) * tau_0 * ( 1.0 - exp(-this->stressGrowthRate * gamma) ) / gamma2 +
                    tau_0 *this->stressGrowthRate *exp(-this->stressGrowthRate *gamma) / gamma;
            mu = mu_0 + tau_0 * ( 1. - exp(-this->stressGrowthRate * gamma) ) / gamma;

            dgde.at(1) = 2.0 * epsd.at(1) / gamma;
            dgde.at(2) = 2.0 * epsd.at(2) / gamma;
            dgde.at(3) = 2.0 * epsd.at(3) / gamma;
            dgde.at(4) = 1.0 * epsd.at(4) / gamma;
            dgde.at(5) = 1.0 * epsd.at(5) / gamma;
            dgde.at(6) = 1.0 * epsd.at(6) / gamma;
#ifdef BINGHAM_MIN_SHEAR_RATE
        }
#endif

        for (int i = 1; i <= 6; i++) {
            answer.at(1,i) = fabs( 2.0 * epsd.at(1) * dmudg * dgde.at(i) );
            answer.at(2,i) = fabs( 2.0 * epsd.at(2) * dmudg * dgde.at(i) );
            answer.at(3,i) = fabs( 2.0 * epsd.at(3) * dmudg * dgde.at(i) );
            answer.at(4,i) = fabs( epsd.at(4) * dmudg * dgde.at(i) );
            answer.at(5,i) = fabs( epsd.at(5) * dmudg * dgde.at(i) );
            answer.at(6,i) = fabs( epsd.at(6) * dmudg * dgde.at(i) );
        }

        answer.at(1,1) += 2.0 * mu;
        answer.at(2,2) += 2.0 * mu;
        answer.at(3,3) += 2.0 * mu;
        answer.at(4,4) += mu;
        answer.at(5,5) += mu;
        answer.at(6,6) += mu;

        return;

    }  else {
        _error("giveDeviatoricStiffnessMatrix: unsupportted material mode");
    }
}


int
BinghamFluidMaterial2 :: checkConsistency()
{
    if ( domain->giveEngngModel()->giveEquationScalingFlag() ) {
        double scale;
        scale = domain->giveEngngModel()->giveVariableScale(VST_Density);
        propertyDictionary->at('d') /= scale;

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
#ifdef BINGHAM_MIN_SHEAR_RATE
    shearRate = max(shearRate, BINGHAM_MIN_SHEAR_RATE);
#endif
    return ( mu_0 + tau_0 * ( 1. - exp(-this->stressGrowthRate * shearRate) ) / shearRate );

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
        _val = (1/6.)*((epsd.at(1)-epsd.at(2))*(epsd.at(1)-epsd.at(2))+ epsd.at(3)*epsd.at(3)/4.0 ) ;
    } else if ( mmode == _2dAxiFlow ) {
        _val = (1/6.)*((epsd.at(1)-epsd.at(2))*(epsd.at(1)-epsd.at(2)) + (epsd.at(2)-epsd.at(3))*(epsd.at(2)-epsd.at(3)) +
                (epsd.at(3)-epsd.at(1))*(epsd.at(3)-epsd.at(1)) + epsd.at(4) * epsd.at(4)/4.0);
    } else if ( mmode == _3dFlow ) {
        _val = 2.0 * (epsd.at(1) * epsd.at(1) + epsd.at(2) * epsd.at(2) + epsd.at(3) * epsd.at(3) )
                    + epsd.at(4) * epsd.at(4) + epsd.at(5) * epsd.at(5) + epsd.at(6) * epsd.at(6);
    } else {
        _error("computeDevStrainMagnitude: unsupported material mode");
    }

    return sqrt(_val);
}

double
BinghamFluidMaterial2 :: computeDevStressMagnitude(MaterialMode mmode, const FloatArray &sigd)
{
    double _val = 0.0;
    if ( mmode == _2dFlow ) {
        _val = (1/6.)*((sigd.at(1)-sigd.at(2))*(sigd.at(1)-sigd.at(2))+ sigd.at(3)*sigd.at(3)) ;
    } else if ( mmode == _2dAxiFlow ) {
        _val = (1/6.)*((sigd.at(1)-sigd.at(2))*(sigd.at(1)-sigd.at(2)) + (sigd.at(2)-sigd.at(3))*(sigd.at(2)-sigd.at(3)) +
                (sigd.at(3)-sigd.at(1))*(sigd.at(3)-sigd.at(1)) + sigd.at(4) * sigd.at(4));

    } else if ( mmode == _3dFlow ) {
        _val = 0.5 * ( sigd.at(1) * sigd.at(1) + sigd.at(2) * sigd.at(2) + sigd.at(3) * sigd.at(3) +
                2.0 * sigd.at(4) * sigd.at(4) + 2.0 * sigd.at(5) * sigd.at(5) + 2.0 * sigd.at(6) * sigd.at(6) );
    } else {
        _error("computeDevStrainMagnitude: unsupported material mode");
    }

    return sqrt(_val);
}

void
BinghamFluidMaterial2 :: computeDeviatoricStrain(FloatArray &answer, const FloatArray &eps, MaterialMode mmode)
{
    if ( mmode ==  _2dFlow ) {
        double ekk=(eps.at(1)+eps.at(2))/3.0;
        //double ekk = 0.0;

        answer.resize(3);
        answer.at(1) = eps.at(1) - ekk;
        answer.at(2) = eps.at(2) - ekk;
        answer.at(3) = eps.at(3);
    } else if ( mmode == _2dAxiFlow ) {
        double ekk=(eps.at(1)+eps.at(2)+eps.at(3))/3.0;
        //double ekk = 0.0;

        answer.resize(4);
        answer.at(1) = eps.at(1) - ekk;
        answer.at(2) = eps.at(2) - ekk;
        answer.at(3) = eps.at(3) - ekk;
        answer.at(4) = eps.at(4);
    } else if ( mmode == _3dFlow ) {
        //double ekk=(eps.at(1)+eps.at(2)+eps.at(3))/3.0;
        double ekk = 0.0;

        answer.resize(6);
        answer.at(1) = eps.at(1) - ekk;
        answer.at(2) = eps.at(2) - ekk;
        answer.at(3) = eps.at(3) - ekk;
        answer.at(4) = eps.at(4);
        answer.at(5) = eps.at(5);
        answer.at(6) = eps.at(6);
    } else {
        _error("computeDeviatoricStrain: unsupported material mode");
    }
}


void
BinghamFluidMaterial2 :: computeDeviatoricStress(FloatArray &answer, const FloatArray &deps,
        double _nu, MaterialMode mmode)
{
    if ( mmode == _2dFlow ) {
        answer.at(1) = 2.0 * _nu * ( deps.at(1) );
        answer.at(2) = 2.0 * _nu * ( deps.at(2) );
        answer.at(3) = deps.at(3) * _nu;
    } else if ( mmode == _2dAxiFlow ) {
        answer.at(1) = 2.0 * _nu * ( deps.at(1) );
        answer.at(2) = 2.0 * _nu * ( deps.at(2) );
        answer.at(3) = 2.0 * _nu * ( deps.at(3) );
        answer.at(4) = deps.at(4) * _nu;
    } else if ( mmode == _3dFlow ) {
        answer.at(1) = 2.0 * _nu * ( deps.at(1) );
        answer.at(2) = 2.0 * _nu * ( deps.at(2) );
        answer.at(3) = 2.0 * _nu * ( deps.at(3) );
        answer.at(4) = deps.at(4) * _nu;
        answer.at(5) = deps.at(5) * _nu;
        answer.at(6) = deps.at(6) * _nu;
    } else {
        _error("computeDeviatoricStrain: unsupported material mode");
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
        _error("BinghamFluidMaterial2Status: unsupported material mode");
    }

    deviatoricStrainVector.resize(_size);
    deviatoricStrainVector.zero();
    deviatoricStressVector.resize(_size);
    deviatoricStressVector.zero();
    temp_deviatoricStrainVector = deviatoricStrainVector;
}

void
BinghamFluidMaterial2Status :: printOutputAt(FILE *File, TimeStep *tNow)
// Prints the strains and stresses on the data file.
{
    int i, n;

    fprintf(File, " strains ");
    n = deviatoricStrainVector.giveSize();
    for ( i = 1; i <= n; i++ ) {
        fprintf( File, " % .4e", deviatoricStrainVector.at(i) );
    }

    fprintf(File, "\n deviatoric stresses");
    n = deviatoricStressVector.giveSize();
    for ( i = 1; i <= n; i++ ) {
        fprintf( File, " % .4e", deviatoricStressVector.at(i) );
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
    deviatoricStrainVector = temp_deviatoricStrainVector;

    // huhu (forcing initial viscosity at the beginning of each time step)
    temp_devStrainMagnitude = 0.0;
    temp_devStressMagnitude = 0.0;
}


void
BinghamFluidMaterial2Status :: initTempStatus()
//
// initialize record at the beginning of new load step
//
{
    FluidDynamicMaterialStatus :: initTempStatus();

    temp_devStrainMagnitude = devStrainMagnitude;
    temp_devStressMagnitude = devStressMagnitude;
    temp_deviatoricStrainVector = deviatoricStrainVector;
}


contextIOResultType
BinghamFluidMaterial2Status :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full ms context (saves state variables, that completely describe
// current state)
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( ( iores = FluidDynamicMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->write(& devStrainMagnitude, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& devStressMagnitude, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = deviatoricStrainVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
BinghamFluidMaterial2Status :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
{
    // FloatArray *s;
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( ( iores = FluidDynamicMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->read(& devStrainMagnitude, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& devStressMagnitude, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = deviatoricStrainVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    return CIO_OK;
}



void
BinghamFluidMaterial2 :: __debug(GaussPoint *gp, TimeStep *atTime)
{
    FloatArray eps_i(3), eps(3), tau(3), tau_p(3), tau_t(3);

    FloatMatrix d;
    int i, nincr = 10000;

    eps_i.at(1) = 0.00005;
    eps_i.at(2) = -0.00002;
    eps_i.at(3) = 0.00002;

    for ( i = 1; i <= nincr; i++ ) {
        eps.add(eps_i);
        computeDeviatoricStressVector(tau, gp, eps, atTime);
        giveDeviatoricStiffnessMatrix(d, TangentStiffness, gp, atTime);
        tau_t.beProductOf(d, eps_i);
        tau_t.add(tau_p);
        //tau.printYourself();
        //tau_t.printYourself();
        //d.printYourself();
        printf( "%e %e %e  %e %e %e %e %e %e\n", eps.at(1), eps.at(2), eps.at(3), tau.at(1), tau.at(2), tau.at(3), tau_t.at(1), tau_t.at(2), tau_t.at(3) );
        tau_p = tau_t;
    }
}

int
BinghamFluidMaterial2 :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    BinghamFluidMaterial2Status *status = ( ( BinghamFluidMaterial2Status * ) this->giveStatus(aGaussPoint) );
    if (type == IST_DeviatoricStressMeasure) {
        answer.setValues(1, status->giveDevStressMagnitude());
        return 1;
    } else {
        return FluidDynamicMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}

InternalStateValueType
BinghamFluidMaterial2 :: giveIPValueType(InternalStateType type)
{
    if (type == IST_DeviatoricStressMeasure ) {
        return ISVT_SCALAR;
    } else {
        return FluidDynamicMaterial :: giveIPValueType(type);
    }
}


int
BinghamFluidMaterial2 :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if (type == IST_DeviatoricStressMeasure ) {
        answer.resize(1); answer.at(1) = 1;
        return 1;
    }  else {
        return FluidDynamicMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}


int
BinghamFluidMaterial2 :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if (type == IST_DeviatoricStressMeasure ) {
        return 1;
    } else {
        return FluidDynamicMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}


} // end namespace oofem
