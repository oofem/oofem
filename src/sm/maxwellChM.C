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

#include "mathfem.h"
#include "maxwellChM.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "gausspnt.h"
#include "structuralcrosssection.h"
#include "timestep.h"
#include "contextioerr.h"

namespace oofem {
MaxwellChainMaterial :: MaxwellChainMaterial(int n, Domain *d) : RheoChainMaterial(n, d)
{}


void
MaxwellChainMaterial :: computeCharCoefficients(FloatArray &answer, GaussPoint *gp,
                                                double atTime)
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
                                            gp, rTimes,
                                            atTime,
                                            atTime);

    // assemble the matrix of the set of linear equations
    // for computing the optimal moduli
    for ( i = 1; i <= this->nUnits; i++ ) {
        taui = this->giveCharTime(i);
        for ( j = 1; j <= this->nUnits; j++ ) {
            tauj = this->giveCharTime(j);
            for ( sum = 0., r = 1; r <= rSize; r++ ) {
                tti = pow( ( atTime + rTimes.at(r) ) / taui, giveCharTimeExponent(i) ) -
                      pow( atTime / taui, giveCharTimeExponent(i) );
                ttj = pow( ( atTime + rTimes.at(r) ) / tauj, giveCharTimeExponent(j) ) -
                      pow( atTime / tauj, giveCharTimeExponent(j) );
                sum += exp(-tti - ttj);
            }

            A.at(i, j) = sum;
        }

        // assemble rhs
        for ( sumRhs = 0., r = 1; r <= rSize; r++ ) {
            tti = pow( ( atTime + rTimes.at(r) ) / taui, giveCharTimeExponent(i) ) -
                  pow( atTime / taui, giveCharTimeExponent(i) );
            sumRhs += exp(-tti) * discreteRelaxFunctionVal.at(r);
        }

        rhs.at(i) = sumRhs;
    }

    // solve the linear system
    A.solveForRhs(rhs, answer);
}



double
MaxwellChainMaterial :: giveEModulus(GaussPoint *gp, TimeStep *atTime)
{
    /*
     * This function returns the incremental modulus for the given time increment.
     * The modulus may also depend on the specimen geometry (gp - dependence).
     *
     * Note: time -1 refers to the previous time.
     */
    int mu;
    double lambdaMu, Emu, deltaYmu;
    double E = 0.0;

    this->updateEparModuli(gp, relMatAge + ( atTime->giveTargetTime() - 0.5 * atTime->giveTimeIncrement() ) / timeFactor);
    for ( mu = 1; mu <= nUnits; mu++ ) {
        deltaYmu = atTime->giveTimeIncrement() / timeFactor / this->giveCharTime(mu);
        if ( deltaYmu <= 0.0 ) {
            deltaYmu = 1.e-3;
        }

        deltaYmu = pow( deltaYmu, this->giveCharTimeExponent(mu) );

        lambdaMu = ( 1.0 - exp(-deltaYmu) ) / deltaYmu;
        Emu      = this->giveEparModulus(mu); // previously updated by updateEparModuli
        E += lambdaMu * Emu;
    }

    return E;
}




void
MaxwellChainMaterial :: giveEigenStrainVector(FloatArray &answer, MatResponseForm form,
                                              GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
//
// computes the strain due to creep at constant stress during the increment
// (in fact, the INCREMENT of creep strain is computed for mode == VM_Incremental)
//
{
    int mu;
    double E;
    double deltaYmu;
    FloatArray *sigmaMu, help, reducedAnswer;
    FloatMatrix B;
    MaxwellChainMaterialStatus *status = ( MaxwellChainMaterialStatus * ) this->giveStatus(gp);

    if ( mode == VM_Incremental ) {
        this->giveUnitComplianceMatrix(B, ReducedForm, gp, atTime);
        reducedAnswer.resize( B.giveNumberOfRows() );

        for ( mu = 1; mu <= nUnits; mu++ ) {
            deltaYmu = atTime->giveTimeIncrement() / timeFactor / this->giveCharTime(mu);
            deltaYmu = pow( deltaYmu, this->giveCharTimeExponent(mu) );
            sigmaMu  = status->giveHiddenVarsVector(mu);
            if ( sigmaMu ) {
                help.beProductOf(B, * sigmaMu); // B can be moved before sum !!!
                help.times( 1.0 - exp(-deltaYmu) );
                reducedAnswer.add(help);
            }
        }

        E = this->giveEModulus(gp, atTime);
        reducedAnswer.times(1.0 / E);

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
MaxwellChainMaterial :: updateYourself(GaussPoint *gp, TimeStep *tNow)
{
    /*
     * Updates hidden variables used to effectively trace the load history
     */


    int mu;
    double deltaYmu, Emu, lambdaMu;
    FloatArray help, *muthHiddenVarsVector, deltaEps0, help1;
    FloatMatrix Binv;
    MaxwellChainMaterialStatus *status =
        ( MaxwellChainMaterialStatus * ) this->giveStatus(gp);

    this->giveUnitStiffnessMatrix(Binv, ReducedForm, gp, tNow);
    help = status->giveTempStrainVector();
    help.subtract( status->giveStrainVector() );

    // Subtract the stress-independent part of strain
    this->computeTrueStressIndependentStrainVector(deltaEps0, gp, tNow, VM_Incremental);
    if ( deltaEps0.giveSize() ) {
        help.subtract(deltaEps0);
    }

    help1.beProductOf(Binv, help);

    this->updateEparModuli(gp, relMatAge + ( tNow->giveTargetTime() - 0.5 * tNow->giveTimeIncrement() ) / timeFactor);

    for ( mu = 1; mu <= nUnits; mu++ ) {
        deltaYmu = tNow->giveTimeIncrement() / timeFactor / this->giveCharTime(mu);
        deltaYmu = pow( deltaYmu, this->giveCharTimeExponent(mu) );

        lambdaMu = ( 1.0 - exp(-deltaYmu) ) / deltaYmu;
        Emu      = this->giveEparModulus(mu);

        muthHiddenVarsVector = status->giveHiddenVarsVector(mu);
        help = help1;
        help.times(lambdaMu * Emu);
        if ( muthHiddenVarsVector ) {
            muthHiddenVarsVector->times( exp(-deltaYmu) );
            muthHiddenVarsVector->add(help);
        } else {
            status->letHiddenVarsVectorBe( mu, (new FloatArray(help)) );
        }
    }

    // now we call MaxwellChainMaterialStatus->updateYourself()
    status->updateYourself(tNow);
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
    RheoChainMaterial :: initializeFrom(ir);
    return IRRT_OK;
}

/****************************************************************************************/

MaxwellChainMaterialStatus :: MaxwellChainMaterialStatus(int n, Domain *d,
                                                         GaussPoint *g, int nunits) :
    RheoChainMaterialStatus(n, d, g, nunits) {}


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
MaxwellChainMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( stream == NULL ) {
        _error("saveContext : can't write into NULL stream");
    }

    if ( ( iores = RheoChainMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
MaxwellChainMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores = RheoChainMaterialStatus :: restoreContext(stream, mode, obj);

    if ( iores != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} // end namespace oofem
