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
#include "rheoChM.h"
#include "material.h"
#include "isolinearelasticmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "structuralcrosssection.h"
#include "contextioerr.h"

namespace oofem {
RheoChainMaterial :: RheoChainMaterial(int n, Domain *d) : StructuralMaterial(n, d),
    EparVal(), charTimes(), discreteTimeScale()
{
    nUnits = 0;
    relMatAge = 0.0;
    linearElasticMaterial = NULL;
    EparValTime = -1.0;
}


RheoChainMaterial :: ~RheoChainMaterial()
{
    if ( linearElasticMaterial ) {
        delete linearElasticMaterial;
    }
}


int
RheoChainMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether the receiver supports the given mode
//
{
    return mode == _3dMat || mode == _PlaneStress ||
        mode == _PlaneStrain || mode == _1dMat ||
        mode == _PlateLayer || mode == _2dBeamLayer;
}


void
RheoChainMaterial :: giveRealStressVector(FloatArray &answer, 
                                          GaussPoint *gp,
                                          const FloatArray &totalStrain,
                                          TimeStep *atTime)
//
// returns the total stress vector (in full or reduced form - form parameter)
// of the receiver according to
// the previous level of stress and the current
// strain increment at time atTime.
//
{
    FloatArray stressIncrement, stressVector, strainIncrement, reducedStrain, eigenStrain;
    FloatMatrix Binv;
    double Emodulus;
    RheoChainMaterialStatus *status = static_cast< RheoChainMaterialStatus * >( this->giveStatus(gp) );

    // initialize the temporary material status and the Gauss point
    this->initTempStatus(gp);
    this->initGpForNewStep(gp);

    // subtract the part of strain which does not depend on the stress increment from the final strain
    // note: it is somewhat confusing to say "stress-dependent" part of strain
    // because here we subtract not only the part of strain that does not depend
    // on stress (i.g., thermal and shrinkage strain) but also the creep strain,
    // which depends on the previous stress history but not on the stress increment
    //
    this->giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain,
                                                atTime, VM_Incremental);

    // subtract the initial strain to get the "net" strain increment,
    // which is related to the stress increment
    strainIncrement.beDifferenceOf(reducedStrain, status->giveStrainVector());

    // get the initial stress, or set it to zero if not available (first step)
    if ( status->giveStressVector().giveSize() ) {
        stressVector = status->giveStressVector();
    } else {
        stressVector.resize( strainIncrement.giveSize() );
        stressVector.zero();
    }

    // evaluate the incremental modulus
    Emodulus = this->giveEModulus(gp, atTime);
    // construct the unit stiffness matrix (depends on Poisson's ratio)
    this->giveUnitStiffnessMatrix(Binv, gp, atTime);
    // multiply the "net" strain increment by unit stiffness and by the incremental modulus
    stressIncrement.beProductOf(Binv, strainIncrement);
    stressIncrement.times(Emodulus);
    // add the stress increment to the inital stress
    stressVector.add(stressIncrement);

    // store the final strain and stress in the status
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(stressVector);

    // update the shrinkage strain if needed
    if ( this->hasIncrementalShrinkageFormulation() ) {
        FloatArray shv;
        this->giveShrinkageStrainVector(shv, gp, atTime, VM_Total);
        status->setShrinkageStrainVector(shv);
    }

    answer = stressVector;
}


void
RheoChainMaterial :: computeDiscreteRelaxationFunction(FloatArray &answer,
                                                       const FloatArray &atTimes,
                                                       double t0, double tr)
{
    int size, nsteps, si;
    double taui, tauk, Jtrt0, totalDeltaSigma;
    double sig0;
    double sum;

    size = atTimes.giveSize();
    nsteps = size;
    FloatArray deltaSigma(size);
    answer.resize(size);
    answer.zero();

    Jtrt0 = this->computeCreepFunction(t0 + atTimes.at(1), t0);
    sig0  = 1. / Jtrt0;
    answer.at(1) = sig0;
    si = 2;

    totalDeltaSigma = 0.;
    for ( int k = si; k <= nsteps; k++ ) {
        sum = 0.;
        for ( int i = si; i <= k - 1; i++ ) {
            if ( i == 1 ) {
                // ??? it seems that tr plays no role because si=2 and the case
                // i=1 or k=1 can never occur
                taui = 0.5 * ( t0 + atTimes.at(i) + tr );
            } else {
                taui = t0 + 0.5 * ( atTimes.at(i) + atTimes.at(i - 1) );
            }

            sum += deltaSigma.at(i) * this->computeCreepFunction(t0 + atTimes.at(k), taui);
        }

        if ( k == 1 ) {
            // ??? it seems that tr plays no role because si=2 and the case
            // i=1 or k=1 can never occur
            tauk = 0.5 * ( t0 + atTimes.at(k) + tr );
        } else {
            tauk = t0 + 0.5 * ( atTimes.at(k) + atTimes.at(k - 1) );
        }

        deltaSigma.at(k) = ( sig0 * ( this->computeCreepFunction(t0 + atTimes.at(k), t0) - Jtrt0 )
                             - sum ) / this->computeCreepFunction(t0 + atTimes.at(k), tauk);

        totalDeltaSigma += deltaSigma.at(k);
        answer.at(k) = sig0 - totalDeltaSigma;
    }
}


void
RheoChainMaterial :: generateLogTimeScale(FloatArray &answer, double from, double to, int nsteps)
{
    answer.resize(nsteps);
    answer.zero();

    double help = log(to/from) / nsteps;
    for ( int i = 1; i <= nsteps; i++ ) {
        answer.at(i) = exp( i * help ) * from;
    }
}


const FloatArray &
RheoChainMaterial :: giveDiscreteTimes()
{
    return discreteTimeScale;
}


void
RheoChainMaterial :: giveUnitStiffnessMatrix(FloatMatrix &answer,
                                             GaussPoint *gp,
                                             TimeStep *tStep)
{
    /*
     * Returns the stiffness matrix of an isotropic linear elastic material
     * with the given Poisson ratio (assumed to be unaffected by creep)
     * and a unit value of Young's modulus.
     *
     * We call crossSection, because we need Binv matrix also for modes
     * defined at the crosssection level.
     * (Hidden stresses are stored in MatStatus level, but they are
     * defined by stressStrain mode in gp, which may be of a type suported
     * at the crosssection level).
     */

    static_cast< StructuralCrossSection * >( gp->giveCrossSection() )
    ->giveCharMaterialStiffnessMatrix(answer, ElasticStiffness, gp, tStep);
}

void
RheoChainMaterial :: giveUnitComplianceMatrix(FloatMatrix &answer,
                                              GaussPoint *gp,
                                              TimeStep *tStep)
/*
 * Returns the compliance matrix of an isotropic linear elastic material
 * with the given Poisson ratio (assumed to be unaffected by creep)
 * and a unit value of Young's modulus.
 */
{
    FloatMatrix tangent;
    static_cast< StructuralCrossSection * >( gp->giveCrossSection() )->
        giveCharMaterialStiffnessMatrix(tangent, ElasticStiffness, gp, tStep);
    answer.beInverseOf(tangent);
}



double
RheoChainMaterial :: giveEparModulus(int iChain)
{
    /* returns the modulus of unit number iChain,
     * previously computed by updateEparModuli() function
     */
    return EparVal.at(iChain);
}


void
RheoChainMaterial :: updateEparModuli(double atTime)
{
    /*
     * Computes moduli of individual units in the chain that provide
     * the best approximation of the relaxation or creep function,
     * depending on whether a Maxwell or Kelvin chain is used.
     *
     * INPUTS:
     *
     * atTime - age of material when load is applied ???
     *
     * DESCRIPTION:
     * We store the computed values because they will be used by other material points in subsequent
     * calculations. Their computation is very costly.
     * For a given active time step, there should be requests only for Epar values at time
     * (currentTime+prevTime)*0.5. The previous values are not needed.
     *
     */
    // compute new values and store them in a temporary array for further use
    if ( fabs(atTime - EparValTime) > TIME_DIFF ) {
        if ( atTime < 0 ) {
            this->computeCharCoefficients(EparVal, 1.e-3);
        } else {
            this->computeCharCoefficients(EparVal, atTime);
        }
        EparValTime = atTime;
    }
}


void
RheoChainMaterial :: computeTrueStressIndependentStrainVector(FloatArray &answer,
                                                              GaussPoint *gp, TimeStep *stepN, ValueModeType mode)
//
// computes the strain due to temperature and shrinkage effects
// (it is called "true" because the "stress-independent strain"
//  will also contain the effect of creep)
//
{
    FloatArray e0;

    // shrinkage strain
    this->giveShrinkageStrainVector(answer, gp, stepN, mode);
    // thermally induced strain
    StructuralMaterial :: computeStressIndependentStrainVector(e0, gp, stepN, mode);
    answer.add(e0);
}


void
RheoChainMaterial :: computeStressIndependentStrainVector(FloatArray &answer,
                                                          GaussPoint *gp, TimeStep *stepN, ValueModeType mode)
//
// computes the strain due to temperature, shrinkage and creep effects
// (it is the strain which would occur at the end of the step if the stress
//  during the step was held constant, so it is not really the part of strain
//  independent of stress, but independent of the stress increment)
//
// takes into account the form of the load vector assumed by engngModel (Incremental or Total Load form)
//
{
    FloatArray e0;

    // strain due to temperature changes and shrinkage
    this->computeTrueStressIndependentStrainVector(answer, gp, stepN, mode);
    // strain due to creep
    this->giveEigenStrainVector(e0, gp, stepN, mode);
    if ( e0.giveSize() ) {
        answer.add(e0);
    }
}


double
RheoChainMaterial :: giveCharTime(int i) const
{
    return charTimes.at(i);
}

void
RheoChainMaterial :: computeCharTimes()
{
    /*
     * This function generates discrete characteristic times
     * according to rules which guarantee a good approximation
     * of the relaxation or creep function by the Dirichlet series
     *
     * Default value of the first relaxation time Tau(1) is chosen to be equal to 0.1 day.
     * The last relaxation time Tau(n) is chosen as 1.0e30
     * (to approximate very long processes).
     * The second largest time Tau(n-1) is chosen as 0.75 tmax
     * where tmax is the lifetime of structure or the end of the time of interest.
     * Times  Tau(2) .. Tau(n-2) are defined by uniform division to n-2 steps in the
     * log scale. It is necessary to check the condition a <= 10, where Tau(k) = a Tau(k-1)
     *
     *
     */

    int size, nsteps;
    double endTime, Taun1, Tau1, help;

    endTime = this->giveEndOfTimeOfInterest() + relMatAge;
    Taun1 = 0.75 * endTime;

    if ( this->begOfTimeOfInterest == -1 ) {
        this->begOfTimeOfInterest = 0.1;         //default value
    }

    Tau1  = begOfTimeOfInterest;

    if ( Tau1 <= 0 ) {
        _error("begOfTimeOfInterest must be a positive number");
    }

    nsteps = ( int ) ( ( log(Taun1) - log(Tau1) ) / log(10.) + 1. );
    if ( nsteps < 8 ) {
        nsteps = 8;
    }

    this->nUnits = size = nsteps + 2;

    this->charTimes.resize(size);

    charTimes.at(1) = Tau1;
    charTimes.at(size) = 1.0e10;
    help = ( log(Taun1) - log(Tau1) ) / ( double ) nsteps;
    for ( int i = 1; i <= nsteps; i++ ) {
        charTimes.at(i + 1) = exp(log(Tau1) + help * i);
    }
}


void
RheoChainMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                   MatResponseMode mode,
                                                   GaussPoint *gp,
                                                   TimeStep *atTime)
{
    //
    // Returns the incremental material stiffness matrix of the receiver
    //
    this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, atTime);
    if ( mode == ElasticStiffness ) return;
    answer.times( this->giveEModulus(gp, atTime) );
}


void
RheoChainMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer,
                                              MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *atTime)
{
    //
    // Returns the incremental material stiffness matrix of the receiver
    //
    this->giveLinearElasticMaterial()->givePlaneStressStiffMtrx(answer, mode, gp, atTime);
    if ( mode == ElasticStiffness ) return;
    answer.times( this->giveEModulus(gp, atTime) );
}

void
RheoChainMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                              MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *atTime)
{
    //
    // Returns the incremental material stiffness matrix of the receiver
    //
    this->giveLinearElasticMaterial()->givePlaneStrainStiffMtrx(answer, mode, gp, atTime);
    if ( mode == ElasticStiffness ) return;
    answer.times( this->giveEModulus(gp, atTime) );
}


void
RheoChainMaterial :: give1dStressStiffMtrx(FloatMatrix &answer,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime)
{
    //
    // Returns the incremental material stiffness matrix of the receiver
    //
    this->giveLinearElasticMaterial()->give1dStressStiffMtrx(answer, mode, gp, atTime);
    if ( mode == ElasticStiffness ) return;
    answer.times( this->giveEModulus(gp, atTime) );
}


MaterialStatus *
RheoChainMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates a new material status corresponding to this class
 */
{
    return new RheoChainMaterialStatus(1, this->giveDomain(), gp, nUnits);
}


IRResultType
RheoChainMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    StructuralMaterial :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, nu, _IFT_RheoChainMaterial_n);
    IR_GIVE_FIELD(ir, relMatAge, _IFT_RheoChainMaterial_relmatage);
    this->begOfTimeOfInterest = -1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, begOfTimeOfInterest, _IFT_RheoChainMaterial_begoftimeofinterest);
    this->endOfTimeOfInterest = -1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, endOfTimeOfInterest, _IFT_RheoChainMaterial_endoftimeofinterest);
    IR_GIVE_FIELD(ir, timeFactor, _IFT_RheoChainMaterial_timefactor); // solution time/timeFactor should give time in days

    // sets up nUnits variable
    this->computeCharTimes();
    // sets up discrete times
    double endTime = this->giveEndOfTimeOfInterest();
    this->generateLogTimeScale(discreteTimeScale, this->begOfTimeOfInterest, endTime, MNC_NPOINTS - 1);
    ///@warning Stiffness is time dependant, so the variable changes with time.
    this->updateEparModuli(0.); // stiffnesses are time independent (evaluated at time t = 0.)

    return IRRT_OK;
}

LinearElasticMaterial *
RheoChainMaterial :: giveLinearElasticMaterial()
{
    if ( linearElasticMaterial == NULL ) {
        linearElasticMaterial = new IsotropicLinearElasticMaterial(this->giveNumber(),
                                                                   this->giveDomain(),
                                                                   1.0, this->nu);
    }

    return linearElasticMaterial;
}

double
RheoChainMaterial :: giveEndOfTimeOfInterest()
{
    if ( this->endOfTimeOfInterest > 0.0 ) {
        return this->endOfTimeOfInterest;
    } else {
        this->endOfTimeOfInterest = this->giveDomain()->giveEngngModel()->giveEndOfTimeOfInterest() / timeFactor;
    }

    return this->endOfTimeOfInterest;
}

contextIOResultType
RheoChainMaterial :: saveIPContext(DataStream *stream, ContextMode mode, GaussPoint* gp)
//
// saves full status for this material, also invokes saving
// for sub-objects of this
{
    contextIOResultType iores;

    if ( stream == NULL ) {
        _error("saveContext : can't write into NULL stream");
    }

    if ( ( iores = Material :: saveIPContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = giveLinearElasticMaterial()->saveIPContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
RheoChainMaterial :: restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp)
//
// reads full status for this material, also invokes reading
// of sub-objects of this
//
{
    contextIOResultType iores;

    if ( ( iores = Material :: restoreIPContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // invoke possible restoring of statuses for submaterials
    if ( ( iores = giveLinearElasticMaterial()->restoreIPContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}



/****************************************************************************************/

RheoChainMaterialStatus :: RheoChainMaterialStatus(int n, Domain *d,
                                                   GaussPoint *g, int nunits) :
    StructuralMaterialStatus(n, d, g), shrinkageStrain() {
    // constructor
    nUnits = nunits;

    hiddenVars = new FloatArray * [ nUnits ];
    for ( int i = 0; i < nUnits; i++ ) {
        hiddenVars [ i ] = NULL;
    }
}


RheoChainMaterialStatus :: ~RheoChainMaterialStatus()
{
    // destructor

    if ( hiddenVars ) {
        for ( int i = 0; i < nUnits; i++ ) {
            delete hiddenVars [ i ];
        }

        delete [] hiddenVars;
    }
}

FloatArray *
RheoChainMaterialStatus :: letHiddenVarsVectorBe(int i, FloatArray *newVector)
{
    /*
     * sets i-th  hidden variables vector to newVector
     * deletes the previous one, if defined
     */
    if ( i > nUnits ) {
        OOFEM_ERROR("RheoChainMaterialStatus :: letHiddenVarsVectorBe - unit number exceeds the specified limit");
    }

    if ( hiddenVars [ i - 1 ] ) {
        delete hiddenVars [ i - 1 ];
    }

    hiddenVars [ i - 1 ] = newVector;

    return newVector;
}



void
RheoChainMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    // printing of hidden variables
    // useful only for debugging
    FloatArray helpVec;

    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "{hidden variables: ");
    for ( int i = 0; i < nUnits; i++ ) {
        StructuralMaterial :: giveFullSymVectorForm(helpVec, * ( hiddenVars [ i ] ), gp->giveMaterialMode());
        fprintf(file, "{ ");
        for ( int j = 1; j <= helpVec.giveSize(); j++ ) {
            fprintf( file, "%f ", helpVec.at(j) );
        }

        fprintf(file, "} ");
    }

    if ( shrinkageStrain.isNotEmpty() ) {
        fprintf(file, "shrinkageStrain: {");
        for ( int j = 1; j <= shrinkageStrain.giveSize(); j++ ) {
            fprintf( file, "%f ", shrinkageStrain.at(j) );
        }

        fprintf(file, "} ");
    }

    fprintf(file, "}\n");
}


void
RheoChainMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
}

void
RheoChainMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
}

contextIOResultType
RheoChainMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
//
{
    contextIOResultType iores;

    if ( stream == NULL ) {
        _error("saveContext : can't write into NULL stream");
    }

    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data
    for ( int i = 0; i < nUnits; i++ ) {
        if ( hiddenVars [ i ] == NULL ) {
            hiddenVars [ i ] = new FloatArray(0);
        }

        if ( ( iores = hiddenVars [ i ]->storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( ( iores = shrinkageStrain.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
RheoChainMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore the state variables from a stream
//
{
    contextIOResultType iores;

    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    for ( int i = 0; i < nUnits; i++ ) {
        if ( hiddenVars [ i ] == NULL ) {
            hiddenVars [ i ] = new FloatArray(0);
        }

        if ( ( iores = hiddenVars [ i ]->restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( ( iores = shrinkageStrain.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} // end namespace oofem
