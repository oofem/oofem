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
#include "kelvinChSolM.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "timestep.h"
#include "contextioerr.h"


namespace oofem {
KelvinChainSolidMaterial :: KelvinChainSolidMaterial(int n, Domain *d) : RheoChainMaterial(n, d)
{ }

double
KelvinChainSolidMaterial :: giveEModulus(GaussPoint *gp, TimeStep *tStep)
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

    if ( this->EparVal.isEmpty() ) {
      this->updateEparModuli(0.); // stiffnesses are time independent (evaluated at time t = 0.)
    }


    for ( mu = 1; mu <= nUnits; mu++ ) {
        lambdaMu = this->computeLambdaMu(gp, tStep, mu);
        Emu = this->giveEparModulus(mu);
        sum += ( 1 - lambdaMu ) / Emu;
    }

    v = this->computeSolidifiedVolume(gp, tStep);
    return sum / v;
}

void
KelvinChainSolidMaterial :: giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode)
//
// computes the strain due to creep at constant stress during the increment
// (in fact, the INCREMENT of creep strain is computed for mode == VM_Incremental)
//
{
    int mu;
    double betaMu;
    double v;
    FloatArray *sigmaVMu = NULL, reducedAnswer, help;
    FloatMatrix C;
    KelvinChainSolidMaterialStatus *status = static_cast< KelvinChainSolidMaterialStatus * >( this->giveStatus(gp) );


    if ( this->EparVal.isEmpty() ) {
      this->updateEparModuli(0.); // stiffnesses are time independent (evaluated at time t = 0.)
    }


    if ( mode == VM_Incremental ) {
        for ( mu = 1; mu <= nUnits; mu++ ) {
            betaMu = this->computeBetaMu(gp, tStep, mu);
            sigmaVMu =  & status->giveHiddenVarsVector(mu); // JB

            if ( sigmaVMu->isNotEmpty() ) {
                help.zero();
                help.add(* sigmaVMu);
                help.times( ( 1.0 - betaMu ) / this->giveEparModulus(mu) );
                reducedAnswer.add(help);
            }
        }

        if ( sigmaVMu->isNotEmpty() ) {
            help = reducedAnswer;
            this->giveUnitComplianceMatrix(C, gp, tStep);
            reducedAnswer.beProductOf(C, help);
            v = this->computeSolidifiedVolume(gp, tStep);
            reducedAnswer.times(1. / v);
        }

        answer = reducedAnswer;
    } else {
        /* error - total mode not implemented yet */
        OOFEM_ERROR("mode is not supported");
    }
}

double
KelvinChainSolidMaterial :: computeBetaMu(GaussPoint *gp, TimeStep *tStep, int Mu)
{
    double betaMu;
    double deltaT;
    double tauMu;

    deltaT = tStep->giveTimeIncrement() / timeFactor;
    tauMu = this->giveCharTime(Mu);

    if ( deltaT / tauMu > 30 ) {
        betaMu = 0;
    } else {
        betaMu = exp(-( deltaT ) / tauMu);
    }

    return betaMu;
}

double
KelvinChainSolidMaterial :: computeLambdaMu(GaussPoint *gp, TimeStep *tStep, int Mu)
{
    double lambdaMu;
    double deltaT;
    double tauMu;

    deltaT = tStep->giveTimeIncrement() / timeFactor;
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
KelvinChainSolidMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    RheoChainMaterial :: giveRealStressVector(answer, gp, reducedStrain, tStep);

    // Computes hidden variables and stores them as temporary
    this->computeHiddenVars(gp, tStep);
}


void
KelvinChainSolidMaterial :: computeHiddenVars(GaussPoint *gp, TimeStep *tStep)
{
    /*
     * Updates hidden variables used to effectively trace the load history
     */

    double betaMu;
    double lambdaMu;

    FloatArray help, SigmaVMu, deltaEps0, deltaSigma;
    FloatMatrix D;
    KelvinChainSolidMaterialStatus *status = static_cast< KelvinChainSolidMaterialStatus * >( this->giveStatus(gp) );


    if ( !this->isActivated(tStep) ) {
      help.resize(StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
      help.zero();
      for ( int mu = 1; mu <= nUnits; mu++ ) {
	status->letTempHiddenVarsVectorBe(mu, help);
      }
      return;
    }
    
    help = status->giveTempStrainVector(); // gives updated strain vector (at the end of time-step)
    help.subtract( status->giveStrainVector() ); // strain increment in current time-step

    // Subtract the stress-independent part of strain
    this->computeStressIndependentStrainVector(deltaEps0, gp, tStep, VM_Incremental);
    if ( deltaEps0.giveSize() ) {
        help.subtract(deltaEps0); // should be equal to zero if there is no stress change during the time-step
    }

    this->giveUnitStiffnessMatrix(D, gp, tStep);

    help.times( this->giveEModulus(gp, tStep) );
    deltaSigma.beProductOf(D, help);

    for ( int mu = 1; mu <= nUnits; mu++ ) {
        betaMu = this->computeBetaMu(gp, tStep, mu);
        lambdaMu = this->computeLambdaMu(gp, tStep, mu);

        help = deltaSigma;
        help.times(lambdaMu);

        SigmaVMu = status->giveHiddenVarsVector(mu);

        if ( SigmaVMu.giveSize() ) {
            SigmaVMu.times(betaMu);
            SigmaVMu.add(help);
            status->letTempHiddenVarsVectorBe(mu, SigmaVMu);
        } else {
            status->letTempHiddenVarsVectorBe(mu, help);
        }
    }
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
    return RheoChainMaterial :: initializeFrom(ir);
}

// useless here
double
KelvinChainSolidMaterial :: computeCreepFunction(double tStep, double ofAge)
{
    OOFEM_ERROR("function has not been yet implemented to KelvinChainSolidMaterialStatus.C");
    return 0.;
}


/****************************************************************************************/

KelvinChainSolidMaterialStatus :: KelvinChainSolidMaterialStatus(int n, Domain *d,
                                                                 GaussPoint *g, int nunits) :
    RheoChainMaterialStatus(n, d, g, nunits) { }

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
KelvinChainSolidMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = RheoChainMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
KelvinChainSolidMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = RheoChainMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} // end namespace oofem
