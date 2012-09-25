/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
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

#include "mathfem.h"
#include "kelvinChM.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "gausspnt.h"
#include "structuralcrosssection.h"
#include "timestep.h"
#include "contextioerr.h"


namespace oofem {
KelvinChainMaterial :: KelvinChainMaterial(int n, Domain *d) : RheoChainMaterial(n, d)
{}

void
KelvinChainMaterial :: computeCharCoefficients(FloatArray &answer, GaussPoint *gp,
                                               double atTime)
{
    /*
     * This function computes the moduli of individual Kelvin units
     * such that the corresponding Dirichlet series gives the best
     * approximation of the actual compliance function J(t,t0) with t0 fixed.
     *
     * The optimal moduli are obtained using the least-square method.
     *
     * INPUTS:
     * atTime = age of material when load is applied
     */

    int i, j, r, rSize;
    double taui, tauj, sum, tti, ttj, sumRhs;
    FloatArray rhs(this->nUnits);
    FloatMatrix A(this->nUnits, this->nUnits);

    const FloatArray &rTimes = this->giveDiscreteTimes();
    rSize = rTimes.giveSize();
    FloatArray discreteComplianceFunctionVal(rSize);

    // compute values of the compliance function at specified times rTimes
    // (can be done directly, since the compliance function is available)

    for ( i = 1; i <= rSize; i++ ) {
        discreteComplianceFunctionVal.at(i) = this->computeCreepFunction(gp, atTime + rTimes.at(i), atTime);
    }

    // assemble the matrix of the set of linear equations
    // for computing the optimal compliances
    // !!! chartime exponents are assumed to be equal to 1 !!!
    for ( i = 1; i <= this->nUnits; i++ ) {
        taui = this->giveCharTime(i);
        for ( j = 1; j <= this->nUnits; j++ ) {
            tauj = this->giveCharTime(j);
            for ( sum = 0., r = 1; r <= rSize; r++ ) {
                tti = rTimes.at(r) / taui;
                ttj = rTimes.at(r) / tauj;
                sum += ( 1. - exp(-tti) ) * ( 1. - exp(-ttj) );
            }

            A.at(i, j) = sum;
        }

        // assemble rhs
        // !!! chartime exponents are assumed to be equal to 1 !!!
        for ( sumRhs = 0., r = 1; r <= rSize; r++ ) {
            tti = rTimes.at(r) / taui;
            sumRhs += ( 1. - exp(-tti) ) * discreteComplianceFunctionVal.at(r);
        }

        rhs.at(i) = sumRhs;
    }

    // solve the linear system
    A.solveForRhs(rhs, answer);

    // convert compliances into moduli
    for ( i = 1; i <= this->nUnits; i++ ) {
        answer.at(i) = 1. / answer.at(i);
    }
}

double
KelvinChainMaterial :: giveEModulus(GaussPoint *gp, TimeStep *atTime)
{
    /*
     * This function returns the incremental modulus for the given time increment.
     * Return value is the incremental E modulus of non-aging Kelvin chain without the first unit (elastic spring)
     * The modulus may also depend on the specimen geometry (gp - dependence).
     *
     * Note: time -1 refers to the previous time.
     */

    int mu;
    double deltaT, tauMu, lambdaMu, Dmu;
    double sum = 0.0; // return value

    if ( EparVal.giveSize() == 0 ) {
        this->updateEparModuli(gp, relMatAge + ( atTime->giveTargetTime() - 0.5 * atTime->giveTimeIncrement() ) / timeFactor);
    }

    deltaT = atTime->giveTimeIncrement() / timeFactor;

    // EparVal values were determined using the least-square method
    for ( mu = 1; mu <= nUnits; mu++ ) {
        tauMu = this->giveCharTime(mu);
        if ( deltaT / tauMu < 1.e-5 ) {
            lambdaMu = 1 - 0.5 * ( deltaT / tauMu ) + 1 / 6 * ( pow(deltaT / tauMu, 2) ) - 1 / 24 * ( pow(deltaT / tauMu, 3) );
        } else if ( deltaT / tauMu > 30 ) {
            lambdaMu = tauMu / deltaT;
        } else {
            lambdaMu = ( 1.0 - exp(-deltaT / tauMu) ) * tauMu / deltaT;
        }

        Dmu = this->giveEparModulus(mu);
        sum += ( 1 - lambdaMu ) / Dmu;
    }

    return sum;
}

void
KelvinChainMaterial :: giveEigenStrainVector(FloatArray &answer, MatResponseForm form,
                                             GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
//
// computes the strain due to creep at constant stress during the increment
// (in fact, the INCREMENT of creep strain is computed for mode == VM_Incremental)
//
{
    int mu;
    double beta;
    FloatArray *gamma, reducedAnswer, help;
    KelvinChainMaterialStatus *status = ( KelvinChainMaterialStatus * ) this->giveStatus(gp);

    // !!! chartime exponents are assumed to be equal to 1 !!!

    if ( mode == VM_Incremental ) {
        for ( mu = 1; mu <= nUnits; mu++ ) {
            if ( ( atTime->giveTimeIncrement() / timeFactor ) / this->giveCharTime(mu) > 30 ) {
                beta = 0;
            } else {
                beta = exp( -( atTime->giveTimeIncrement() / timeFactor ) / ( this->giveCharTime(mu) ) );
            }

            gamma =  status->giveHiddenVarsVector(mu);
            if ( gamma ) {
                help.zero();
                help.add(* gamma);
                help.times(1.0 - beta);
                reducedAnswer.add(help);
            }
        }

        if ( form == ReducedForm ) {
            answer =  reducedAnswer;
            return;
        }

        // expand the strain to full form if requested
        ( ( StructuralCrossSection * ) gp->giveCrossSection() )->
        giveFullCharacteristicVector(answer, gp, reducedAnswer);
    } else {
        /* error - total mode not implemented yet */
        _error("giveEigenStrainVector - mode is not supported");
    }
}


void
KelvinChainMaterial :: updateYourself(GaussPoint *gp, TimeStep *tNow)
{
    /*
     * Updates hidden variables used to effectively trace the load history
     */

    // !!! chartime exponents are assumed to be equal to 1 !!!
    double betaMu;
    double lambdaMu;
    double deltaT;
    double tauMu;

    FloatArray help, *muthHiddenVarsVector, deltaEps0, help1;
    KelvinChainMaterialStatus *status = ( KelvinChainMaterialStatus * ) this->giveStatus(gp);

    help = status->giveTempStrainVector(); // gives updated strain vector (at the end of time-step)
    help.subtract( status->giveStrainVector() ); // strain increment in current time-step

    // Subtract the stress-independent part of strain
    this->computeStressIndependentStrainVector(deltaEps0, gp, tNow, VM_Incremental);
    if ( deltaEps0.giveSize() ) {
        help.subtract(deltaEps0); // should be equal to zero if there is no stress change during the time-step
    }

    help.times( this->giveEModulus(gp, tNow) ); // = delta_sigma
    help1 = help;

    deltaT = tNow->giveTimeIncrement() / timeFactor;

    for ( int mu = 1; mu <= nUnits; mu++ ) {
        help = help1;
        tauMu = ( this->giveCharTime(mu) );

        muthHiddenVarsVector = status->giveHiddenVarsVector(mu); //gamma_mu

        if ( deltaT / tauMu < 1.e-5 ) {
            betaMu = exp(-( deltaT ) / tauMu);
            lambdaMu = 1 - 0.5 * ( deltaT / tauMu ) + 1 / 6 * ( pow(deltaT / tauMu, 2) ) - 1 / 24 * ( pow(deltaT / tauMu, 3) );
        } else if ( deltaT / tauMu > 30 )       {
            betaMu = 0;
            lambdaMu = tauMu / deltaT;
        } else   {
            betaMu = exp(-( deltaT ) / tauMu);
            lambdaMu = ( 1.0 - betaMu ) * tauMu / deltaT;
        }

        help.times( lambdaMu / ( this->giveEparModulus(mu) ) );

        if ( muthHiddenVarsVector ) {
            muthHiddenVarsVector->times(betaMu);
            muthHiddenVarsVector->add(help);
        } else {
            status->letHiddenVarsVectorBe( mu, (new FloatArray(help)) );
        }
    }

    // now we call KelvinChainMaterialStatus->updateYourself()
    status->updateYourself(tNow);
}



MaterialStatus *
KelvinChainMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates a new material status corresponding to this class
 */
{
    return new KelvinChainMaterialStatus(1, this->giveDomain(), gp, nUnits);
}

IRResultType
KelvinChainMaterial :: initializeFrom(InputRecord *ir)
{
    RheoChainMaterial :: initializeFrom(ir);
    this->giveDiscreteTimes(); // Makes sure the new discrete times are evaluated.
    return IRRT_OK;
}

/****************************************************************************************/

KelvinChainMaterialStatus :: KelvinChainMaterialStatus(int n, Domain *d,
                                                       GaussPoint *g, int nunits) :
    RheoChainMaterialStatus(n, d, g, nunits) {}

void
KelvinChainMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    RheoChainMaterialStatus :: printOutputAt(file, tStep);
}


void
KelvinChainMaterialStatus :: updateYourself(TimeStep *tStep)
{
    RheoChainMaterialStatus :: updateYourself(tStep);
}

void
KelvinChainMaterialStatus :: initTempStatus()
{
    RheoChainMaterialStatus :: initTempStatus();
}

contextIOResultType
KelvinChainMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = RheoChainMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
KelvinChainMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = RheoChainMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} // end namespace oofem
