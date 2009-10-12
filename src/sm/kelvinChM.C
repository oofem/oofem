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
 *               Copyright (C) 1993 - 2009   Borek Patzak
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

// file: KelvinChM.C

#ifndef __MAKEDEPEND
#include <math.h>
#endif
#include "mathfem.h"
#include "kelvinChM.h"
#include "material.h"
#include "isolinearelasticmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "matstatus.h"
#include "gausspnt.h"
#include "structuralcrosssection.h"
#include "timestep.h"
#include "contextioerr.h"


KelvinChainMaterial :: KelvinChainMaterial(int n, Domain *d) : RheoChainMaterial(n, d)
{
}

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
    FloatArray rhs(this->nUnits), discreteComplianceFunctionVal;
    FloatMatrix A(this->nUnits, this->nUnits);

    const FloatArray &rTimes = this->giveDiscreteTimes();
    rSize = rTimes.giveSize();

    // compute values of the compliance function at specified times rTimes
    // (can be done directly, since the compliance function is available)
    for ( i = 1; i<= rSize; j++)
      discreteComplianceFunctionVal.at(i) = this->computeCreepFunction(gp, atTime + rTimes.at(i), atTime);

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
                sum += (1.-exp(-tti)) * (1.-exp(-ttj));
            }

            A.at(i, j) = sum;
        }

        // assemble rhs
	// !!! chartime exponents are assumed to be equal to 1 !!!
        for ( sumRhs = 0., r = 1; r <= rSize; r++ ) {
            tti = rTimes.at(r) / taui;
            sumRhs += (1.-exp(-tti)) * discreteComplianceFunctionVal.at(r);
        }

        rhs.at(i) = sumRhs;
    }

    // solve the linear system
    A.solveForRhs(rhs, answer);

    // convert compliances into moduli
    for ( i = 1; i <= this->nUnits; i++ )
      answer.at(i) = 1./answer.at(i);

    return;
}



double
KelvinChainMaterial :: giveEModulus(GaussPoint *gp, TimeStep *atTime)
{
    /*
     * This function returns the incremental modulus for the given time increment.
     * The modulus may also depend on the specimen geometry (gp - dependence).
     *
     * It is stored as "Einc" for further expected requests from other gaussPoints that correspond to the same material.
     *
     * Note: time -1 refers to the previous time.    
     */
    int mu;
    double lambdaMu, Emu, deltaYmu;
    double E = 0.0;

    // !!! this should be replaced by the expression valid for the Kelvin chain !!!
    // !!! chartime exponents are assumed to be equal to 1 !!!
    this->updateEparModuli(gp, relMatAge + ( atTime->giveTime() - 0.5 * atTime->giveTimeIncrement() ) / timeFactor);
    for ( mu = 1; mu <= nUnits; mu++ ) {
        deltaYmu = atTime->giveTimeIncrement() / timeFactor / this->giveCharTime(mu);
        if ( deltaYmu <= 0.0 ) {
            deltaYmu = 1.e-3;
        }

        deltaYmu = __OOFEM_POW( deltaYmu, this->giveCharTimeExponent(mu) );

        lambdaMu = ( 1.0 - exp(-deltaYmu) ) / deltaYmu;
        Emu      = this->giveEparModulus(mu); // previously updated by updateEparModuli

        E += lambdaMu * Emu;
    }

    Einc = E;
    return Einc;
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
    double E;
    double deltaYmu;
    FloatArray *sigmaMu, help, reducedAnswer;
    FloatMatrix B;
    KelvinChainMaterialStatus *status = ( KelvinChainMaterialStatus * ) this->giveStatus(gp);

    // !!! this should be replaced by the expression valid for the Kelvin chain !!!
    // !!! chartime exponents are assumed to be equal to 1 !!!
    if ( mode == VM_Incremental ) {
        this->giveUnitComplianceMatrix(B, ReducedForm, gp, atTime);
        reducedAnswer.resize( B.giveNumberOfRows() );

        for ( mu = 1; mu <= nUnits; mu++ ) {
            deltaYmu = atTime->giveTimeIncrement() / timeFactor / this->giveCharTime(mu);
            deltaYmu = __OOFEM_POW( deltaYmu, this->giveCharTimeExponent(mu) );
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

    return;
}


void
KelvinChainMaterial :: updateYourself(GaussPoint *gp, TimeStep *tNow)
{
    /*
     * Updates hidden variables used to effectively trace the load history
     */
    // !!! this should be replaced by the expression valid for the Kelvin chain !!!
    // !!! chartime exponents are assumed to be equal to 1 !!!
    int mu;
    double deltaYmu, Emu, lambdaMu;
    FloatArray help, *muthHiddenVarsVector, deltaEps0, help1;
    FloatMatrix Binv;
    KelvinChainMaterialStatus *status =
        ( KelvinChainMaterialStatus * ) this->giveStatus(gp);

    this->giveUnitStiffnessMatrix(Binv, ReducedForm, gp, tNow);
    help = status->giveTempStrainVector();
    help.substract( status->giveStrainVector() );

    // Subtract the stress-independent part of strain
    this->computeTrueStressIndependentStrainVector(deltaEps0, gp, tNow, VM_Incremental);
    if ( deltaEps0.giveSize() ) {
        help.substract(deltaEps0);
    }

    help1.beProductOf(Binv, help);

    this->updateEparModuli(gp, relMatAge + ( tNow->giveTime() - 0.5 * tNow->giveTimeIncrement() ) / timeFactor);

    for ( mu = 1; mu <= nUnits; mu++ ) {
        deltaYmu = tNow->giveTimeIncrement() / timeFactor / this->giveCharTime(mu);
        deltaYmu = __OOFEM_POW( deltaYmu, this->giveCharTimeExponent(mu) );

        lambdaMu = ( 1.0 - exp(-deltaYmu) ) / deltaYmu;
        Emu      = this->giveEparModulus(mu);

        muthHiddenVarsVector = status->giveHiddenVarsVector(mu);
        help = help1;
        help.times(lambdaMu * Emu);
        if ( muthHiddenVarsVector ) {
            muthHiddenVarsVector->times( exp(-deltaYmu) );
            muthHiddenVarsVector->add(& help);
        } else {
            status->letHiddenVarsVectorBe( mu, help.GiveCopy() );
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
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    RheoChainMaterial :: initializeFrom(ir);

    return IRRT_OK;
}

contextIOResultType
KelvinChainMaterial :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = RheoChainMaterial :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
KelvinChainMaterial :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = RheoChainMaterial :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}



/****************************************************************************************/

KelvinChainMaterialStatus :: KelvinChainMaterialStatus(int n, Domain *d,
                                                         GaussPoint *g, int nunits) :
  RheoChainMaterialStatus(n, d, g, nunits) {
}

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
    int i;
    contextIOResultType iores;

    if ( ( iores = RheoChainMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    return CIO_OK;
}

