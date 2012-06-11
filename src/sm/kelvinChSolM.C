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
#include "kelvinChSolM.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "gausspnt.h"
#include "structuralcrosssection.h"
#include "timestep.h"
#include "contextioerr.h"


namespace oofem {
KelvinChainSolidMaterial :: KelvinChainSolidMaterial(int n, Domain *d) : RheoChainMaterial(n, d)
{}

double
KelvinChainSolidMaterial :: giveEModulus(GaussPoint *gp, TimeStep *atTime)
{
    /*
     * This function returns the incremental modulus for the given time increment.
     * Return value is the incremental E modulus of non-aging Kelvin chain without the first unit (elastic spring)
     * The modulus may also depend on the specimen geometry (gp - dependence).
     *
     * It is stored as "Einc" for further expected requests from other gaussPoints that correspond to the same material.
     *
     * Note: time -1 refers to the previous time.
     */

    int mu;
    double v;
    double lambdaMu, Emu;
    double sum = 0.0;

    if ( EparVal.giveSize() == 0 ) {
        this->updateEparModuli(gp, 0.); // stiffnesses are time independent (evaluated at time t = 0.)
    }

    for ( mu = 1; mu <= nUnits; mu++ ) {
        lambdaMu = this->computeLambdaMu(gp, atTime, mu);
        Emu = this->giveEparModulus(mu);
        sum += ( 1 - lambdaMu ) / Emu;
    }

    v = this->computeSolidifiedVolume(gp, atTime);
    return sum / v;
}

void
KelvinChainSolidMaterial :: giveEigenStrainVector(FloatArray &answer, MatResponseForm form,
                                                  GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
//
// computes the strain due to creep at constant stress during the increment
// (in fact, the INCREMENT of creep strain is computed for mode == VM_Incremental)
//
{
    int mu;
    double betaMu;
    double v;
    FloatArray *sigmaVMu, reducedAnswer, help;
    FloatMatrix C;
    KelvinChainSolidMaterialStatus *status = ( KelvinChainSolidMaterialStatus * ) this->giveStatus(gp);

    if ( mode == VM_Incremental ) {
        for ( mu = 1; mu <= nUnits; mu++ ) {
            betaMu = this->computeBetaMu(gp, atTime, mu);
            sigmaVMu =  status->giveHiddenVarsVector(mu);

            if ( sigmaVMu ) {
                help.zero();
                help.add(* sigmaVMu);
                help.times( ( 1.0 - betaMu ) / this->giveEparModulus(mu) );
                reducedAnswer.add(help);
            }
        }

        if ( sigmaVMu ) {
            help = reducedAnswer;
            this->giveUnitComplianceMatrix(C, ReducedForm, gp, atTime);
            reducedAnswer.beProductOf(C, help);
            v = this->computeSolidifiedVolume(gp, atTime);
            reducedAnswer.times(1. / v);
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

double
KelvinChainSolidMaterial :: computeBetaMu(GaussPoint *gp, TimeStep *atTime, double Mu)
{
    double betaMu;
    double deltaT;
    double tauMu;

    deltaT = atTime->giveTimeIncrement() / timeFactor;
    tauMu = this->giveCharTime(Mu);

    if ( deltaT / tauMu > 30 ) {
        betaMu = 0;
    } else {
        betaMu = exp(-( deltaT ) / tauMu);
    }

    return betaMu;
}

double
KelvinChainSolidMaterial :: computeLambdaMu(GaussPoint *gp, TimeStep *atTime, double Mu)
{
    double lambdaMu;
    double deltaT;
    double tauMu;

    deltaT = atTime->giveTimeIncrement() / timeFactor;
    tauMu = this->giveCharTime(Mu);

    if ( deltaT / tauMu < 1.e-5 ) {
        lambdaMu = 1 - 0.5 * ( deltaT / tauMu ) + 1 / 6 * ( pow(deltaT / tauMu, 2) ) - 1 / 24 * ( pow(deltaT / tauMu, 3) );
    } else if ( deltaT / tauMu > 30 ) {
        lambdaMu = tauMu / deltaT;
    } else {
        lambdaMu = ( 1.0 -  exp(-( deltaT ) / tauMu) ) * tauMu / deltaT;
    }

    return lambdaMu;
}

void
KelvinChainSolidMaterial :: updateYourself(GaussPoint *gp, TimeStep *tNow)
{
    /*
     * Updates hidden variables used to effectively trace the load history
     */

    double betaMu;
    double lambdaMu;

    FloatArray help, *SigmaVMu, deltaEps0, deltaSigma;
    FloatMatrix D;
    KelvinChainSolidMaterialStatus *status = ( KelvinChainSolidMaterialStatus * ) this->giveStatus(gp);

    help = status->giveTempStrainVector(); // gives updated strain vector (at the end of time-step)
    help.subtract( status->giveStrainVector() ); // strain increment in current time-step

    // Subtract the stress-independent part of strain
    this->computeStressIndependentStrainVector(deltaEps0, gp, tNow, VM_Incremental);
    if ( deltaEps0.giveSize() ) {
        help.subtract(deltaEps0); // should be equal to zero if there is no stress change during the time-step
    }


    this->giveUnitStiffnessMatrix(D, ReducedForm, gp, tNow);


    help.times( this->giveEModulus(gp, tNow) );
    deltaSigma.beProductOf(D, help);

    for ( int mu = 1; mu <= nUnits; mu++ ) {
        betaMu = this->computeBetaMu(gp, tNow, mu);
        lambdaMu = this->computeLambdaMu(gp, tNow, mu);

        help = deltaSigma;
        help.times(lambdaMu);

        SigmaVMu = status->giveHiddenVarsVector(mu);

        if ( SigmaVMu ) {
            SigmaVMu->times(betaMu);
            SigmaVMu->add(help);
        } else {
            status->letHiddenVarsVectorBe( mu, ( new FloatArray(help) ) );
        }
    }

    status->updateYourself(tNow);
}


MaterialStatus *
KelvinChainSolidMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates a new material status corresponding to this class
 */
{
    return new KelvinChainSolidMaterialStatus(1, this->giveDomain(), gp, nUnits);
}

IRResultType
KelvinChainSolidMaterial :: initializeFrom(InputRecord *ir)
{
    RheoChainMaterial :: initializeFrom(ir);

    return IRRT_OK;
}

// useless here
double
KelvinChainSolidMaterial :: computeCreepFunction(GaussPoint *gp, double atTime, double ofAge)
{
    _error("computeCreepFunction: function has not been yet implemented to KelvinChainSolidMaterialStatus.C");
    return 0.;
}


/****************************************************************************************/

KelvinChainSolidMaterialStatus :: KelvinChainSolidMaterialStatus(int n, Domain *d,
                                                                 GaussPoint *g, int nunits) :
    RheoChainMaterialStatus(n, d, g, nunits) {}

void
KelvinChainSolidMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    RheoChainMaterialStatus :: printOutputAt(file, tStep);
}


void
KelvinChainSolidMaterialStatus :: updateYourself(TimeStep *tStep)
{
    RheoChainMaterialStatus :: updateYourself(tStep);
}

void
KelvinChainSolidMaterialStatus :: initTempStatus()
{
    RheoChainMaterialStatus :: initTempStatus();
}

contextIOResultType
KelvinChainSolidMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = RheoChainMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
KelvinChainSolidMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = RheoChainMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} // end namespace oofem
