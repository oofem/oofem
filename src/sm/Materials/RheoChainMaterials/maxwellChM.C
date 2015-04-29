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

#include "mathfem.h"
#include "maxwellChM.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "timestep.h"
#include "contextioerr.h"

namespace oofem {
MaxwellChainMaterial :: MaxwellChainMaterial(int n, Domain *d) : RheoChainMaterial(n, d)
{ }


void
MaxwellChainMaterial :: computeCharCoefficients(FloatArray &answer, double tStep)
{
    int i, j, r, rSize;
    double taui, tauj, sum, tti, ttj, sumRhs;
    FloatArray rhs(this->nUnits), discreteRelaxFunctionVal;
    FloatMatrix A(this->nUnits, this->nUnits);

    const FloatArray &rTimes = this->giveDiscreteTimes();
    rSize = rTimes.giveSize();

    // compute discrete values of the relaxation function at times rTimes
    // from the creep function (by numerically solving integral equations)
    //
    // a direct call to the relaxation function could be used, if available
    this->computeDiscreteRelaxationFunction(discreteRelaxFunctionVal,
                                            rTimes,
                                            tStep,
                                            tStep);

    // assemble the matrix of the set of linear equations
    // for computing the optimal moduli
    for ( i = 1; i <= this->nUnits; i++ ) {
        taui = this->giveCharTime(i);
        for ( j = 1; j <= this->nUnits; j++ ) {
            tauj = this->giveCharTime(j);
            for ( sum = 0., r = 1; r <= rSize; r++ ) {
                tti = pow( ( tStep + rTimes.at(r) ) / taui, giveCharTimeExponent(i) ) -
                pow( tStep / taui, giveCharTimeExponent(i) );
                ttj = pow( ( tStep + rTimes.at(r) ) / tauj, giveCharTimeExponent(j) ) -
                pow( tStep / tauj, giveCharTimeExponent(j) );
                sum += exp(-tti - ttj);
            }

            A.at(i, j) = sum;
        }

        // assemble rhs
        for ( sumRhs = 0., r = 1; r <= rSize; r++ ) {
            tti = pow( ( tStep + rTimes.at(r) ) / taui, giveCharTimeExponent(i) ) -
            pow( tStep / taui, giveCharTimeExponent(i) );
            sumRhs += exp(-tti) * discreteRelaxFunctionVal.at(r);
        }

        rhs.at(i) = sumRhs;
    }

    // solve the linear system
    A.solveForRhs(rhs, answer);
}



double
MaxwellChainMaterial :: giveEModulus(GaussPoint *gp, TimeStep *tStep)
{
    /*
     * This function returns the incremental modulus for the given time increment.
     * The modulus may also depend on the specimen geometry (gp - dependence).
     *
     * Note: time -1 refers to the previous time.
     */
    double lambdaMu, Emu, deltaYmu;
    double E = 0.0;

    ///@warning THREAD UNSAFE!
    this->updateEparModuli(relMatAge + ( tStep->giveTargetTime() - 0.5 * tStep->giveTimeIncrement() ) / timeFactor);

    for ( int mu = 1; mu <= nUnits; mu++ ) {
        deltaYmu = tStep->giveTimeIncrement() / timeFactor / this->giveCharTime(mu);
        if ( deltaYmu <= 0.0 ) {
            deltaYmu = 1.e-3;
        }

        deltaYmu = pow( deltaYmu, this->giveCharTimeExponent(mu) );

        lambdaMu = ( 1.0 - exp(-deltaYmu) ) / deltaYmu;
        Emu = this->giveEparModulus(mu); // previously updated by updateEparModuli
        E += lambdaMu * Emu;
    }

    return E;
}




void
MaxwellChainMaterial :: giveEigenStrainVector(FloatArray &answer,
                                              GaussPoint *gp, TimeStep *tStep, ValueModeType mode)
//
// computes the strain due to creep at constant stress during the increment
// (in fact, the INCREMENT of creep strain is computed for mode == VM_Incremental)
//
{
    int mu;
    double E;
    double deltaYmu;
    FloatArray help, reducedAnswer;
    //FloatArray *sigmaMu;
    FloatArray sigmaMu;
    FloatMatrix B;
    MaxwellChainMaterialStatus *status = static_cast< MaxwellChainMaterialStatus * >( this->giveStatus(gp) );

    if ( mode == VM_Incremental ) {
        this->giveUnitComplianceMatrix(B, gp, tStep);
        reducedAnswer.resize( B.giveNumberOfRows() );
        reducedAnswer.zero();

        for ( mu = 1; mu <= nUnits; mu++ ) {
            deltaYmu = tStep->giveTimeIncrement() / timeFactor / this->giveCharTime(mu);
            deltaYmu = pow( deltaYmu, this->giveCharTimeExponent(mu) );

            sigmaMu  = status->giveHiddenVarsVector(mu); // JB

            if ( sigmaMu.giveSize() ) {
                help.beProductOf(B, sigmaMu); // B can be moved before sum !!!
                help.times( 1.0 - exp(-deltaYmu) );
                reducedAnswer.add(help);
            }
        }

        E = this->giveEModulus(gp, tStep);
        reducedAnswer.times(1.0 / E);

        answer =  reducedAnswer;
    } else {
        /* error - total mode not implemented yet */
        OOFEM_ERROR("mode is not supported");
    }
}


void
MaxwellChainMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    RheoChainMaterial :: giveRealStressVector(answer, gp, reducedStrain, tStep);

    // Computes hidden variables and stores them as temporary
    this->computeHiddenVars(gp, tStep);
}


void
MaxwellChainMaterial :: computeHiddenVars(GaussPoint *gp, TimeStep *tStep)
{
    /*
     * Updates hidden variables used to effectively trace the load history
     */

    double deltaYmu, Emu, lambdaMu;
    FloatArray help, deltaEps0, help1;
    FloatArray muthHiddenVarsVector;

    FloatMatrix Binv;
    MaxwellChainMaterialStatus *status =
        static_cast< MaxwellChainMaterialStatus * >( this->giveStatus(gp) );

    this->giveUnitStiffnessMatrix(Binv, gp, tStep);
    help = status->giveTempStrainVector();
    help.subtract( status->giveStrainVector() );

    // Subtract the stress-independent part of strain
    this->computeTrueStressIndependentStrainVector(deltaEps0, gp, tStep, VM_Incremental);
    if ( deltaEps0.giveSize() ) {
        help.subtract(deltaEps0);
    }

    help1.beProductOf(Binv, help);

    ///@warning THREAD UNSAFE!
    this->updateEparModuli(relMatAge + ( tStep->giveTargetTime() - 0.5 * tStep->giveTimeIncrement() ) / timeFactor);

    for ( int mu = 1; mu <= nUnits; mu++ ) {
        deltaYmu = tStep->giveTimeIncrement() / timeFactor / this->giveCharTime(mu);
        deltaYmu = pow( deltaYmu, this->giveCharTimeExponent(mu) );

        lambdaMu = ( 1.0 - exp(-deltaYmu) ) / deltaYmu;
        Emu = this->giveEparModulus(mu);

        muthHiddenVarsVector = status->giveHiddenVarsVector(mu);
        help = help1;
        help.times(lambdaMu * Emu);
        if ( muthHiddenVarsVector.giveSize() ) {
            muthHiddenVarsVector.times( exp(-deltaYmu) );
            muthHiddenVarsVector.add(help);
            status->letTempHiddenVarsVectorBe(mu, muthHiddenVarsVector);
        } else {
            status->letTempHiddenVarsVectorBe(mu, help);
        }
    }
}


MaterialStatus *
MaxwellChainMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates a new material status corresponding to this class
 */
{
    return new MaxwellChainMaterialStatus(1, this->giveDomain(), gp, nUnits);
}


IRResultType
MaxwellChainMaterial :: initializeFrom(InputRecord *ir)
{
    return RheoChainMaterial :: initializeFrom(ir);
}

/****************************************************************************************/

MaxwellChainMaterialStatus :: MaxwellChainMaterialStatus(int n, Domain *d,
                                                         GaussPoint *g, int nunits) :
    RheoChainMaterialStatus(n, d, g, nunits) { }


void
MaxwellChainMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    RheoChainMaterialStatus :: printOutputAt(file, tStep);
}


void
MaxwellChainMaterialStatus :: updateYourself(TimeStep *tStep)
{
    RheoChainMaterialStatus :: updateYourself(tStep);
}

void
MaxwellChainMaterialStatus :: initTempStatus()
{
    RheoChainMaterialStatus :: initTempStatus();
}

contextIOResultType
MaxwellChainMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = RheoChainMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
MaxwellChainMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores = RheoChainMaterialStatus :: restoreContext(stream, mode, obj);

    if ( iores != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} // end namespace oofem
