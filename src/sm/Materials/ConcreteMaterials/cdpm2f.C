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
 *               Copyright (C) 1993 - 2023   Borek Patzak
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

#include "cdpm2f.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "intarray.h"
#include "datastream.h"
#include "contextioerr.h"
#include "timestep.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "mathfem.h"
#include "classfactory.h"
#include "matconst.h"
#include <limits>

namespace oofem {
    REGISTER_Material(CDPM2F);

    CDPM2FStatus::CDPM2FStatus(GaussPoint *gp) :
        ConcreteDPM2Status(gp), RandomMaterialStatusExtensionInterface()
    {}


    Interface *
    CDPM2FStatus::giveInterface(InterfaceType type)
    {
        if ( type == RandomMaterialStatusExtensionInterfaceType ) {
            return static_cast< RandomMaterialStatusExtensionInterface * >(this);
        } else {
            return nullptr;
        }
    }


    //   ********************************
    //   *** CLASS CDPM2F ***
    //   ********************************

#define IDM_ITERATION_LIMIT 1.e-8

    CDPM2F::CDPM2F(int n, Domain *d) :
        ConcreteDPM2(n, d), RandomMaterialExtensionInterface()
    {}


    void CDPM2F::initializeFrom(const std::shared_ptr<InputRecord> &ir)
    {
        // call the corresponding service for the linear elastic material
        ConcreteDPM2::initializeFrom(ir);
        RandomMaterialExtensionInterface::initializeFrom(ir);


        //default not available
        IR_GIVE_FIELD(ir, lf, _IFT_CDPM2F_Lf);

        //Default not available
        IR_GIVE_FIELD(ir, vf0, _IFT_CDPM2F_Vf0);

        //Default not available
        IR_GIVE_FIELD(ir, df, _IFT_CDPM2F_Df);

        //Default not available
        IR_GIVE_FIELD(ir, ef, _IFT_CDPM2F_Ef);

        //Default 1 MPa. Must be entered because of the units.
        IR_GIVE_FIELD(ir, tau0, _IFT_CDPM2F_Tau0);

        //Default 0.015
        this->beta = 0.015;
        IR_GIVE_OPTIONAL_FIELD(ir, beta, _IFT_CDPM2F_Beta);

        this->f = 0.8;
        IR_GIVE_OPTIONAL_FIELD(ir, f, _IFT_CDPM2F_f);

        //default not available.  Theroetical value very small
        IR_GIVE_FIELD(ir, sm, _IFT_CDPM2F_Sm);

        this->alpha = 1.;
        IR_GIVE_OPTIONAL_FIELD(ir, alpha, _IFT_CDPM2F_Alpha);


        this->xi = 10.;
        IR_GIVE_OPTIONAL_FIELD(ir, this->xi, _IFT_CDPM2F_xi);

        //Precalculate parameters which do not depend on alpha
        this->g = 2. * ( 1. + exp(M_PI * this->f / 2.) ) / ( 4. + pow(this->f, 2.) );

        this->c = this->beta * this->lf / ( 2. * this->df );

        this->deltaUl = this->lf / 2.;

        this->inputParameters = this->computeFibreParameters(this->alpha);

        if ( alpha < inputParameters.alphaMin ) {
            OOFEM_ERROR("Wrong input in CDPM2F. alpha=%e should be larger than alphamin %e\n", alpha, inputParameters.alphaMin);
        }
    }


    CDPM2F::FibreParameters
    CDPM2F::computeFibreParameters(double alpha) const
    {
        FibreParameters fp;
        fp.alpha = alpha;

        //Introduce reduction in vf due to dispersion
        fp.vf = ( 1. + log(alpha) ) * this->vf0;

        fp.eta = this->ef * fp.vf / ( this->eM * ( 1. - fp.vf ) );

        fp.s0 = 1. / 2. * this->g * this->tau0 * fp.vf * ( 1. + fp.eta ) * this->lf / this->df;

        fp.omega = sqrt(4. * ( 1. + fp.eta ) * this->beta * this->tau0 / this->ef);

        fp.k = fp.omega * this->lf / ( 2. * this->df );

        fp.lambda = cosh(fp.k) - 1.;

        fp.deltaStar = ( 2. * this->df ) / this->beta * fp.lambda;

        if ( this->c <= 6 * fp.lambda + 2 ) {  //softening starts at end of debonding
            fp.deltaCu = this->lf * fp.lambda / this->c;
        } else { //softening starts later
            fp.deltaCu = ( 0.5 * this->lf * ( this->c - 2. ) / ( 3. * this->c ) );
        }

        //softening starts at end of debonding
        if ( this->c <= 6 * fp.lambda + 2 ) {
            fp.stressCu = fp.s0 * ( 1. + 2. * fp.lambda ) * pow( ( 2. * fp.lambda / this->c - 1. ), 2. );
        } else {
            fp.stressCu = fp.s0 * 4. * pow( ( c + 1. ), 3. ) / 27. / pow(c, 2.);
        }

        //Calculate alphamin to check alpha input.
        double ftTemp = this->ft * ( 1. - this->yieldTolDamage );

        double z = ftTemp / ( 0.5 * g * tau0 * lf / df * ( ( 1. + beta * fp.deltaCu / df ) * ( pow( ( 1. - 2. * fp.deltaCu / lf ), 2.) ) ) );

        double vfm = ( -( this->eM + this->eM * z ) + sqrt(pow( ( this->eM + this->eM * z ), 2 ) + 4. * ( this->ef - this->eM ) * this->eM * z) ) / ( 2. * ( this->ef - this->eM ) );

        fp.alphaMin = exp( ( vfm - this->vf0 ) / this->vf0 );

        return fp;
    }


    CDPM2F::FibreParameters
    CDPM2F::giveFibreParameters(GaussPoint *gp) const
    {
        double randomFactor;
        if ( !RandomMaterialExtensionInterface::give(alphaDispersion_ID, gp, randomFactor) ) {
            return this->inputParameters;
        }

        //The random field value is a factor on the mean alpha given in the input.
        //Random alpha is limited to [alphaMin, 1]. Since alphaMin depends weakly on alpha, the lower limit is found by fixed point iteration.
        double randomAlpha = min(randomFactor * this->alpha, 1.);
        FibreParameters fp = this->computeFibreParameters(randomAlpha);
        for ( int i = 0; i < 20 && fp.alpha < fp.alphaMin * ( 1. - 1.e-10 ); i++ ) {
            fp = this->computeFibreParameters(fp.alphaMin);
        }

        return fp;
    }


    std::unique_ptr< MaterialStatus >
    CDPM2F::CreateStatus(GaussPoint *gp) const
    {
        return std::make_unique< CDPM2FStatus >(gp);
    }


    MaterialStatus *
    CDPM2F::giveStatus(GaussPoint *gp) const
    {
        if ( !gp->hasMaterialStatus() ) {
            gp->setMaterialStatus(this->CreateStatus(gp) );
            this->_generateStatusVariables(gp);
        }
        return static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    }


    double CDPM2F::computeFibreStress(double delta, const FibreParameters &fp) const
    {
        double fibreStress = 0.;


        if ( delta >= 0 && delta <= fp.deltaStar ) {
            fibreStress = ( 2. / fp.k * ( ( 1. - acosh(1. + fp.lambda * delta / fp.deltaStar) / fp.k ) * sqrt(pow( ( 1. + fp.lambda * delta / fp.deltaStar ), 2.) - 1.) + ( fp.lambda * delta ) / ( fp.k * fp.deltaStar ) ) + this->ap * delta / fp.deltaStar ) * fp.s0;
        } else if ( delta > fp.deltaStar && delta <= this->deltaUl ) {
            fibreStress = ( 1. + beta * delta / df ) * ( pow( ( 1. - 2. * delta / lf ), 2.) ) * fp.s0;

            //Check that fibre stress is not negative. Should not be the case
            if ( fibreStress <= 0. ) {
                OOFEM_WARNING("Fibre stress negative in CDPM2::computeFibreStress\n");
                fibreStress = 0.;
            }
        } else if ( delta > this->deltaUl || delta < 0 ) {
            fibreStress = 0.;
        }

        return fibreStress;
    }


    double CDPM2F::computeMatrixStress(double delta, const FibreParameters &fp) const
    {
        double matrixStress = 0.;
        double ftTemp = this->ft * ( 1. - this->yieldTolDamage );

        if ( this->softeningType == 2 ) {
            //Chao, should this be vf or vf0 here?
            matrixStress = ( 1 - fp.vf ) * ftTemp * exp(-delta / this->wf);
        } else {
            OOFEM_ERROR("concrete softening must be exponential (stype = 2)");
        }

        return matrixStress;
    }

    double CDPM2F::computeCrackOpening(double crackingStrain, const double le, const FibreParameters &fp) const
    {
        double delta = 0., gammaR0 = 0., stressRatio = 0., dStressRatioDDelta = 0.;
        int nite = 0;


        double residual = 0., dResidualDDelta = 0.;


        double ftTemp = this->ft * ( 1. - this->yieldTolDamage );


        double gammaCu = ( 1. - fp.alpha ) * ftTemp / this->eM * this->sm / ( fp.deltaCu * ( 1. - fp.alphaMin ) ) + ( ( fp.alpha - fp.alphaMin ) / ( 1. - fp.alphaMin ) );


        double eCu = fp.deltaCu * gammaCu / this->sm;
        double eUl = this->lf / 2. / le;

        double deltaCuUnloading = fp.deltaCu * ( gammaCu * le - sm ) / ( le - sm );

        if ( crackingStrain >= 0 && crackingStrain <= eCu ) { //pre-preak
            delta = fp.deltaCu * ( 1. - exp( -crackingStrain / ( this->xi * eCu ) ) ) / ( 1. - exp(-eCu  / ( this->xi * eCu ) ) ); //sigmoid delta relation
        } else if ( crackingStrain > eCu && crackingStrain <= eUl ) {
            //initial guess of delta
            delta = fp.deltaCu;

            //Two cases: 1) Unloading occurs in the element. 2) No unloading occurs.
            if ( le > sm/gammaCu ) {
                //Case 1: Unloading in element occurs.

                //Apply Newton method to solve third order equation.

                delta = fp.deltaCu;

                do{
                    if ( nite == 1000 ) {
                        OOFEM_ERROR("Method to compute crack opening in CDPM2F did not converge.")
                    }


                    stressRatio = fp.s0 * ( 1. + this->beta * delta / this->df ) * pow(1. - 2. * delta / this->lf, 2.) / fp.stressCu;

                    dStressRatioDDelta = fp.s0 * this->beta / this->df * pow(1. - 2. * delta / this->lf, 2.) / fp.stressCu -
                                         fp.s0 * ( 1. + this->beta * delta / this->df ) * 2. * ( 1. - 2. * delta / this->lf ) / fp.stressCu * 2. / this->lf;


                    residual = 1. / le * ( delta + ( le / this->sm - 1. ) * deltaCuUnloading * stressRatio ) - crackingStrain;

                    dResidualDDelta = 1. / le + 1. / le * ( le / this->sm - 1. ) * deltaCuUnloading * dStressRatioDDelta;

                    delta -= residual / dResidualDDelta;

                    nite++;
                }while ( fabs(residual) / eCu > 1.e-6 );
            } else {
                //Element so small so that no unloading occurs in element.
                //Are we underestimating fracture energy with this approach. If yes, it must be negligible because so much energy is dissipated in pre-preak.

                gammaR0 = gammaCu * le / sm;

                delta = 1. / ( 4. * ( gammaR0 - 1. ) ) *
                    ( lf * gammaR0 - 2. * fp.deltaCu - sqrt(pow(lf, 2.) * pow(gammaR0, 2.) - 4. * fp.deltaCu * ( lf * gammaR0 - fp.deltaCu ) +
                     8. * ( 1. - gammaR0 ) * crackingStrain * le * ( lf - 2. * fp.deltaCu ) ) );
            }
        } else   {
            delta = le * crackingStrain;
        }

        return delta;
    }


    double CDPM2F::computeStressResidual(double equivStrain, double damage, double kappaOne, double kappaTwo, double le, const FibreParameters &fp) const
    {
        //calculate cracking strain
        double crackingStrain  = kappaOne + damage * kappaTwo;

        double delta = computeCrackOpening(crackingStrain, le, fp);

        double fibreStress = computeFibreStress(delta, fp);

        double matrixStress = computeMatrixStress(delta, fp);

        double residual = ( 1. - damage ) * this->eM * equivStrain - fibreStress - matrixStress;

        return residual;
    }


    double
    CDPM2F::computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double damageOld, double rateFactor, GaussPoint *gp) const
    {
        const FibreParameters fp = this->giveFibreParameters(gp);

        // So that damage does not turn out to be negative if function is entered for equivstrains smaller than e0.
        double ftTemp = this->ft * ( 1. - yieldTolDamage );
        double gammaCu = ( 1. - fp.alpha ) * ftTemp / this->eM * this->sm / ( fp.deltaCu * ( 1. - fp.alphaMin ) ) + ( ( fp.alpha - fp.alphaMin ) / ( 1. - fp.alphaMin ) );

        //Check first if element length is small enough. This needs to be done here because le is caluclated from the principal directions at the onset of cracking.
        if ( le > ( 7. * c - 2. ) * this->sm / ( gammaCu * ( 3. * c - 6. ) ) ) {
            OOFEM_WARNING("element size should be less than %e. Your element size is le = %e\n Local snapback could occur.\n", ( 7. * c - 2. ) * this->sm / ( gammaCu * ( 3. * c - 6. ) ), le);
        }



        /*The main idea of the function is to determine the damage variable with the bisection method, which provides an equivalent 1D stress which is equal to the sum of fibre and concrete stress. The input to the function is the cracking strain. Both fibre and concrete stress are computed from the crack opening. Therefore, we need to calculate the crack opening from the cracking strain.
         * This is carried out in three steps. For a given damage variable:
         * 1) The crack opening is determined from the cracking strain. There are two regions which are separated, namely distributed and localised cracking.
         * 2) Use this crack opening to calculate the fibre stress. Since these expressions are highly nonlinear and implicit, we need to use a bisection method for solving for the damage variable.
         * 3) Compute the concrete stress.
         * 4) Check stress balance. If out of balance go to 1) with a different damage variable.
         */

        double dDamage, residual, residualMid, damageMid, damage = 0.;
        double damageOne = 0.;
        double damageTwo = 1.;


        int nite = 0;
        if ( equivStrain > e0 * ( 1. - yieldTolDamage ) ) { // Check if damage should start
            nite     = 0;

            residual = computeStressResidual(equivStrain, damageOne, kappaOne, kappaTwo, le, fp);
            residualMid = computeStressResidual(equivStrain, damageTwo, kappaOne, kappaTwo, le, fp);

            if ( ( residual * residualMid > 0 ) ) {
                /* //Bisection method cannot find a solution. This can be the case, because the fibre law is too steep at tha start. If this is the reason, then damage is kept constant at this stage. */
                return damageOld;
            }

            damage = residual < 0.0 ? ( dDamage = damageTwo - damageOne, damageOne ) : ( dDamage = damageOne - damageTwo, damageTwo );


            do {
                nite++;
                if ( nite == 100 ) {
                    OOFEM_ERROR("In computeDamageTension: bisection method not converged. residual = %e, residualMid= %e, damage=%e", residual, residualMid, damage);
                }

                damageMid = damage + ( dDamage *= 0.5 );
                residualMid = computeStressResidual(equivStrain, damageMid, kappaOne, kappaTwo, le, fp);

                if ( residualMid <= 0.0 ) {
                    damage = damageMid;
                }
            } while ( fabs(dDamage) > 1.e-10 && residualMid != 0.0 );
        }

        if ( damage > 1. ) {
            damage = 1.;
        }

        if ( damage < 0. || damage < damageOld ) {
            damage = damageOld;
        }

        return damage;
    }
}
