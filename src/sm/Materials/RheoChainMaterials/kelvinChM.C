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

FloatArray
KelvinChainMaterial :: computeCharCoefficients(double tPrime, GaussPoint *gp, TimeStep *tStep) const
{
    /*
     * This function computes the moduli of individual Kelvin units
     * such that the corresponding Dirichlet series gives the best
     * approximation of the actual compliance function J(t,t0) with t0 fixed.
     *
     * The optimal moduli are obtained using the least-square method.
     *
     * INPUTS:
     * tPrime = age of material when load is applied
     */

    FloatArray rhs(this->nUnits);
    FloatMatrix A(this->nUnits, this->nUnits);

    const FloatArray &rTimes = this->giveDiscreteTimes();
    int rSize = rTimes.giveSize();
    FloatArray discreteComplianceFunctionVal(rSize);

    // compute values of the compliance function at specified times rTimes
    // (can be done directly, since the compliance function is available)

    for ( int i = 1; i <= rSize; i++ ) {
        discreteComplianceFunctionVal.at(i) = this->computeCreepFunction(tPrime + rTimes.at(i), tPrime, gp, tStep);
    }

    // assemble the matrix of the set of linear equations
    // for computing the optimal compliances
    // !!! chartime exponents are assumed to be equal to 1 !!!
    for ( int i = 1; i <= this->nUnits; i++ ) {
        double taui = this->giveCharTime(i);
        for ( int j = 1; j <= this->nUnits; j++ ) {
            double tauj = this->giveCharTime(j);
            double sum = 0.;
            for ( int r = 1; r <= rSize; r++ ) {
                double tti = rTimes.at(r) / taui;
                double ttj = rTimes.at(r) / tauj;
                sum += ( 1. - exp(-tti) ) * ( 1. - exp(-ttj) );
            }

            A.at(i, j) = sum;
        }

        // assemble rhs
        // !!! chartime exponents are assumed to be equal to 1 !!!
        double sumRhs = 0.;
        for ( int r = 1; r <= rSize; r++ ) {
            double tti = rTimes.at(r) / taui;
            sumRhs += ( 1. - exp(-tti) ) * discreteComplianceFunctionVal.at(r);
        }

        rhs.at(i) = sumRhs;
    }

    // solve the linear system
    FloatArray answer;
    A.solveForRhs(rhs, answer);

    // convert compliances into moduli
    for ( int i = 1; i <= this->nUnits; i++ ) {
        answer.at(i) = 1. / answer.at(i);
    }
    return answer;
}

double
KelvinChainMaterial :: giveEModulus(GaussPoint *gp, TimeStep *tStep) const
{
    /*
     * This function returns the incremental modulus for the given time increment.
     * Return value is the incremental E modulus of non-aging Kelvin chain without the first unit (elastic spring)
     * The modulus may also depend on the specimen geometry (gp - dependence).
     *
     * Note: time -1 refers to the previous time.
     */

    double sum = 0.0; // return value

    // the viscoelastic material does not exist yet
    if  ( ! Material :: isActivated( tStep ) ) {
      OOFEM_ERROR("Attempted to evaluate E modulus at time lower than casting time");
    }

    ///@warning THREAD UNSAFE!
    double tPrime = this->relMatAge - this->castingTime + ( tStep->giveTargetTime() - 0.5 * tStep->giveTimeIncrement() );
    this->updateEparModuli(tPrime, gp, tStep);

    double deltaT = tStep->giveTimeIncrement();

    // EparVal values were determined using the least-square method
    for ( int mu = 1; mu <= nUnits; mu++ ) {
        double tauMu = this->giveCharTime(mu);
        double lambdaMu;
        if ( deltaT / tauMu < 1.e-5 ) {
            lambdaMu = 1 - 0.5 * ( deltaT / tauMu ) + 1 / 6 * ( pow(deltaT / tauMu, 2) ) - 1 / 24 * ( pow(deltaT / tauMu, 3) );
        } else if ( deltaT / tauMu > 30 ) {
            lambdaMu = tauMu / deltaT;
        } else {
            lambdaMu = ( 1.0 - exp(-deltaT / tauMu) ) * tauMu / deltaT;
        }

        double Dmu = this->giveEparModulus(mu);
        sum += ( 1 - lambdaMu ) / Dmu;
    }

    //    return sum;
    // changed formulation to return stiffness instead of compliance
    return 1. / sum;

}

void
KelvinChainMaterial :: giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const
//
// computes the strain due to creep at constant stress during the increment
// (in fact, the INCREMENT of creep strain is computed for mode == VM_Incremental)
//
{
    auto status = static_cast< KelvinChainMaterialStatus * >( this->giveStatus(gp) );

    // !!! chartime exponents are assumed to be equal to 1 !!!

    if ( ! Material :: isActivated( tStep ) ) {
        OOFEM_ERROR("Attempted to evaluate creep strain for time lower than casting time");
    }

    if ( mode == VM_Incremental ) {
        FloatArray reducedAnswer;

        for ( int mu = 1; mu <= nUnits; mu++ ) {
            double beta;
            if ( tStep->giveTimeIncrement() / this->giveCharTime(mu) > 30 ) {
                beta = 0;
            } else {
                beta = exp( - tStep->giveTimeIncrement() / this->giveCharTime(mu) );
            }

            FloatArray *gamma = & status->giveHiddenVarsVector(mu); // JB
            if ( gamma ) {
                reducedAnswer.add(1.0 - beta, * gamma);
            }
        }

        answer = reducedAnswer;
    } else {
        /* error - total mode not implemented yet */
        OOFEM_ERROR("mode is not supported");
    }
}


void
KelvinChainMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) const
{
    RheoChainMaterial :: giveRealStressVector(answer, gp, reducedStrain, tStep);

    this->computeHiddenVars(gp, tStep);
}

void
KelvinChainMaterial :: computeHiddenVars(GaussPoint *gp, TimeStep *tStep) const
{
    /*
     * Updates hidden variables used to effectively trace the load history
     */

    // !!! chartime exponents are assumed to be equal to 1 !!!

    FloatArray help, delta_sigma;
    KelvinChainMaterialStatus *status = static_cast< KelvinChainMaterialStatus * >( this->giveStatus(gp) );

    // goes there if the viscoelastic material does not exist yet
    if (  ! Material :: isActivated( tStep ) )  {
        help.resize(StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        help.zero();
        for ( int mu = 1; mu <= nUnits; mu++ ) {
            status->letTempHiddenVarsVectorBe(mu, help);
        }
        return;
    }

    delta_sigma = status->giveTempStrainVector(); // gives updated strain vector (at the end of time-step)
    delta_sigma.subtract( status->giveStrainVector() ); // strain increment in current time-step

    // Subtract the stress-independent part of strain
    auto deltaEps0 = this->computeStressIndependentStrainVector(gp, tStep, VM_Incremental);
    if ( deltaEps0.giveSize() ) {
        delta_sigma.subtract(deltaEps0); // should be equal to zero if there is no stress change during the time-step
    }

    // no need to worry about "zero-stiffness" for time < castingTime - this is done above
    delta_sigma.times( this->giveEModulus(gp, tStep) ); // = delta_sigma

    double deltaT = tStep->giveTimeIncrement();

    for ( int mu = 1; mu <= nUnits; mu++ ) {
        double betaMu;
        double lambdaMu;
        help = delta_sigma;
        double tauMu = this->giveCharTime(mu);

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

        FloatArray muthHiddenVarsVector = status->giveHiddenVarsVector(mu); //gamma_mu
        if ( muthHiddenVarsVector.giveSize() ) {
            muthHiddenVarsVector.times(betaMu);
            muthHiddenVarsVector.add(help);
            status->letTempHiddenVarsVectorBe(mu, muthHiddenVarsVector);
        } else {
            status->letTempHiddenVarsVectorBe(mu, help);
        }
    }
}


std::unique_ptr<MaterialStatus> 
KelvinChainMaterial :: CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<KelvinChainMaterialStatus>(gp, nUnits);
}

void
KelvinChainMaterial :: initializeFrom(InputRecord &ir)
{
    RheoChainMaterial :: initializeFrom(ir);
    this->giveDiscreteTimes(); // Makes sure the new discrete times are evaluated.
}

/****************************************************************************************/

KelvinChainMaterialStatus :: KelvinChainMaterialStatus(GaussPoint *g, int nunits) :
    RheoChainMaterialStatus(g, nunits) { }

void
KelvinChainMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
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

void
KelvinChainMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    RheoChainMaterialStatus :: saveContext(stream, mode);
}

void
KelvinChainMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    RheoChainMaterialStatus :: restoreContext(stream, mode);
}
} // end namespace oofem
