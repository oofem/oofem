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

#include "mathfem.h"
#include "kelvinChM.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "timestep.h"
#include "contextioerr.h"


namespace oofem {
KelvinChainMaterial :: KelvinChainMaterial(int n, Domain *d) : RheoChainMaterial(n, d)
{ }

void
KelvinChainMaterial :: computeCharCoefficients(FloatArray &answer, double tStep)
{
    /*
     * This function computes the moduli of individual Kelvin units
     * such that the corresponding Dirichlet series gives the best
     * approximation of the actual compliance function J(t,t0) with t0 fixed.
     *
     * The optimal moduli are obtained using the least-square method.
     *
     * INPUTS:
     * tStep = age of material when load is applied
     */

    int rSize;
    double taui, tauj, tti, ttj;
    FloatArray rhs(this->nUnits);
    FloatMatrix A(this->nUnits, this->nUnits);

    const FloatArray &rTimes = this->giveDiscreteTimes();
    rSize = rTimes.giveSize();
    FloatArray discreteComplianceFunctionVal(rSize);

    // compute values of the compliance function at specified times rTimes
    // (can be done directly, since the compliance function is available)

    for ( int i = 1; i <= rSize; i++ ) {
        discreteComplianceFunctionVal.at(i) = this->computeCreepFunction(tStep + rTimes.at(i), tStep);
    }

    // assemble the matrix of the set of linear equations
    // for computing the optimal compliances
    // !!! chartime exponents are assumed to be equal to 1 !!!
    for ( int i = 1; i <= this->nUnits; i++ ) {
        taui = this->giveCharTime(i);
        for ( int j = 1; j <= this->nUnits; j++ ) {
            tauj = this->giveCharTime(j);
            double sum = 0.;
            for ( int r = 1; r <= rSize; r++ ) {
                tti = rTimes.at(r) / taui;
                ttj = rTimes.at(r) / tauj;
                sum += ( 1. - exp(-tti) ) * ( 1. - exp(-ttj) );
            }

            A.at(i, j) = sum;
        }

        // assemble rhs
        // !!! chartime exponents are assumed to be equal to 1 !!!
        double sumRhs = 0.;
        for ( int r = 1; r <= rSize; r++ ) {
            tti = rTimes.at(r) / taui;
            sumRhs += ( 1. - exp(-tti) ) * discreteComplianceFunctionVal.at(r);
        }

        rhs.at(i) = sumRhs;
    }

    // solve the linear system
    A.solveForRhs(rhs, answer);

    // convert compliances into moduli
    for ( int i = 1; i <= this->nUnits; i++ ) {
        answer.at(i) = 1. / answer.at(i);
    }
}

double
KelvinChainMaterial :: giveEModulus(GaussPoint *gp, TimeStep *tStep)
{
    /*
     * This function returns the incremental modulus for the given time increment.
     * Return value is the incremental E modulus of non-aging Kelvin chain without the first unit (elastic spring)
     * The modulus may also depend on the specimen geometry (gp - dependence).
     *
     * Note: time -1 refers to the previous time.
     */

    double deltaT, tauMu, lambdaMu, Dmu;
    double sum = 0.0; // return value

    ///@warning THREAD UNSAFE!
    this->updateEparModuli(relMatAge + ( tStep->giveTargetTime() - 0.5 * tStep->giveTimeIncrement() ) / timeFactor);

    deltaT = tStep->giveTimeIncrement() / timeFactor;

    // EparVal values were determined using the least-square method
    for ( int mu = 1; mu <= nUnits; mu++ ) {
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
KelvinChainMaterial :: giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode)
//
// computes the strain due to creep at constant stress during the increment
// (in fact, the INCREMENT of creep strain is computed for mode == VM_Incremental)
//
{
    double beta;
    FloatArray *gamma, reducedAnswer, help;
    KelvinChainMaterialStatus *status = static_cast< KelvinChainMaterialStatus * >( this->giveStatus(gp) );

    // !!! chartime exponents are assumed to be equal to 1 !!!

    if ( mode == VM_Incremental ) {
        for ( int mu = 1; mu <= nUnits; mu++ ) {
            if ( ( tStep->giveTimeIncrement() / timeFactor ) / this->giveCharTime(mu) > 30 ) {
                beta = 0;
            } else {
                beta = exp( -( tStep->giveTimeIncrement() / timeFactor ) / ( this->giveCharTime(mu) ) );
            }

            gamma = & status->giveHiddenVarsVector(mu); // JB
            if ( gamma ) {
                help.zero();
                help.add(* gamma);
                help.times(1.0 - beta);
                reducedAnswer.add(help);
            }
        }

        answer = reducedAnswer;
    } else {
        /* error - total mode not implemented yet */
        OOFEM_ERROR("mode is not supported");
    }
}


void
KelvinChainMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    RheoChainMaterial :: giveRealStressVector(answer, gp, reducedStrain, tStep);

    this->computeHiddenVars(gp, tStep);
}

void
KelvinChainMaterial :: computeHiddenVars(GaussPoint *gp, TimeStep *tStep)
{
    /*
     * Updates hidden variables used to effectively trace the load history
     */

    // !!! chartime exponents are assumed to be equal to 1 !!!
    double betaMu;
    double lambdaMu;
    double deltaT;
    double tauMu;

    FloatArray help, deltaEps0, delta_sigma;
    //FloatArray *muthHiddenVarsVector;
    FloatArray muthHiddenVarsVector;
    KelvinChainMaterialStatus *status = static_cast< KelvinChainMaterialStatus * >( this->giveStatus(gp) );

    delta_sigma = status->giveTempStrainVector(); // gives updated strain vector (at the end of time-step)
    delta_sigma.subtract( status->giveStrainVector() ); // strain increment in current time-step

    // Subtract the stress-independent part of strain
    this->computeStressIndependentStrainVector(deltaEps0, gp, tStep, VM_Incremental);
    if ( deltaEps0.giveSize() ) {
        delta_sigma.subtract(deltaEps0); // should be equal to zero if there is no stress change during the time-step
    }

    delta_sigma.times( this->giveEModulus(gp, tStep) ); // = delta_sigma

    deltaT = tStep->giveTimeIncrement() / timeFactor;

    for ( int mu = 1; mu <= nUnits; mu++ ) {
        help = delta_sigma;
        tauMu = ( this->giveCharTime(mu) );

        muthHiddenVarsVector = status->giveHiddenVarsVector(mu); //gamma_mu

        if ( deltaT / tauMu < 1.e-5 ) {
            betaMu = exp(-( deltaT ) / tauMu);
            lambdaMu = 1 - 0.5 * ( deltaT / tauMu ) + 1 / 6 * ( pow(deltaT / tauMu, 2) ) - 1 / 24 * ( pow(deltaT / tauMu, 3) );
        } else if ( deltaT / tauMu > 30 ) {
            betaMu = 0;
            lambdaMu = tauMu / deltaT;
        } else {
            betaMu = exp(-( deltaT ) / tauMu);
            lambdaMu = ( 1.0 - betaMu ) * tauMu / deltaT;
        }

        help.times( lambdaMu / ( this->giveEparModulus(mu) ) );

        if ( muthHiddenVarsVector.giveSize() ) {
            muthHiddenVarsVector.times(betaMu);
            muthHiddenVarsVector.add(help);
            status->letTempHiddenVarsVectorBe(mu, muthHiddenVarsVector);
        } else {
            status->letTempHiddenVarsVectorBe(mu, help);
        }
    }
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
    IRResultType result = RheoChainMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    this->giveDiscreteTimes(); // Makes sure the new discrete times are evaluated.
    return IRRT_OK;
}

/****************************************************************************************/

KelvinChainMaterialStatus :: KelvinChainMaterialStatus(int n, Domain *d,
                                                       GaussPoint *g, int nunits) :
    RheoChainMaterialStatus(n, d, g, nunits) { }

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
KelvinChainMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = RheoChainMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
KelvinChainMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = RheoChainMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} // end namespace oofem
