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
 *               Copyright (C) 1993 - 2015   Borek Patzak
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
#include "rheoChM.h"
#include "material.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "sm/Materials/LatticeMaterials/latticelinearelastic.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "engngm.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "contextioerr.h"

namespace oofem {
RheoChainMaterial :: RheoChainMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{}


RheoChainMaterial :: ~RheoChainMaterial()
{
    if ( linearElasticMaterial ) {
        delete linearElasticMaterial;
    }
}


bool
RheoChainMaterial :: hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _3dMat || mode == _PlaneStress ||
           mode == _PlaneStrain || mode == _1dMat ||
           mode == _PlateLayer || mode == _2dBeamLayer ||
           mode == _1dLattice || mode == _2dLattice ||
      mode == _3dLattice;
}


void
RheoChainMaterial :: giveRealStressVector(FloatArray &answer,
                                          GaussPoint *gp,
                                          const FloatArray &totalStrain,
                                          TimeStep *tStep)
//
// returns the total stress vector (in full or reduced form - form parameter)
// of the receiver according to
// the previous level of stress and the current
// strain increment at time tStep.
//
{
    FloatArray stressIncrement, stressVector, strainIncrement, reducedStrain;
    FloatMatrix Binv;
    double Emodulus;
    RheoChainMaterialStatus *status = static_cast< RheoChainMaterialStatus * >( this->giveStatus(gp) );

    // Initialize the temporary material status and the Gauss point
    // Do not initialize if the the actual time is the same as before
    // this has to be introduced when the time-stepping is variable and
    // some time steps can be run several times
    // It is vital to erase some temporary variables if the time is changed and
    // the time step is still the same
    // In order to guarantee a reasonable performance the temporary status is not
    // initialized in subsequent iterations of the same time step

    if ( status->giveCurrentTime() != tStep->giveTargetTime() ) {
        status->setCurrentTime(tStep->giveTargetTime() );
        this->initTempStatus(gp);
    }

    if ( !this->isActivated(tStep) ) {
        FloatArray zeros;
        zeros.resize(StructuralMaterial :: giveSizeOfVoigtSymVector(gp->giveMaterialMode() ) );
        zeros.zero();
        status->letTempStrainVectorBe(zeros);
        status->letTempStressVectorBe(zeros);
        answer = zeros;
        return;
    }

    // subtract the part of strain which does not depend on the stress increment from the final strain
    // note: it is somewhat confusing to say "stress-dependent" part of strain
    // because here we subtract not only the part of strain that does not depend
    // on stress (e.g., thermal and shrinkage strain) but also the creep strain,
    // which depends on the previous stress history but not on the stress increment
    //
    this->giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain,
                                                tStep, VM_Incremental);
    // subtract the initial strain to get the "net" strain increment,
    // which is related to the stress increment
    strainIncrement.beDifferenceOf(reducedStrain, status->giveStrainVector() );

    // get the initial stress, or set it to zero if not available (first step)
    //if ( status->giveStressVector().giveSize() ) {
    //  stressVector = status->giveStressVector();
    if ( status->giveViscoelasticStressVector().giveSize() ) {
        stressVector = status->giveViscoelasticStressVector();
    } else {
        stressVector.resize(strainIncrement.giveSize() );
        stressVector.zero();
    }


    if ( ( !Material :: isActivated(tStep) ) && ( preCastingTimeMat > 0 ) ) {
        StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(preCastingTimeMat) );
        sMat->giveStiffnessMatrix(Binv, ElasticStiffness, gp, tStep);
    } else {
        // evaluate the incremental modulus
        Emodulus = this->giveEModulus(gp, tStep);
        // construct the unit stiffness matrix (depends on Poisson's ratio)
        this->giveUnitStiffnessMatrix(Binv, gp, tStep);
        // multiply the "net" strain increment by unit stiffness and by the incremental modulus
        Binv.times(Emodulus);
    }

    stressIncrement.beProductOf(Binv, strainIncrement);
    //    stressIncrement.times(Emodulus);

    // add the stress increment to the inital stress
    stressVector.add(stressIncrement);

    // store the final strain and stress in the status
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(stressVector);

    if ( Material :: isActivated(tStep) ) {
        // update the shrinkage strain if needed
        if ( this->hasIncrementalShrinkageFormulation() ) {
            FloatArray shv;
            this->giveShrinkageStrainVector(shv, gp, tStep, VM_Total);
            status->setShrinkageStrainVector(shv);
        }
    }

    answer = stressVector;
}


FloatArrayF< 6 >
RheoChainMaterial :: giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const
//
// returns a FloatArray(6) of initial strain vector
// eps_0 = {exx_0, eyy_0, ezz_0, gyz_0, gxz_0, gxy_0}^T
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    MaterialMode mMode =  gp->giveMaterialMode();
    if ( mMode == _1dLattice || mMode ==  _2dLattice || mMode ==  _3dLattice ) {
        return {
                   this->talpha, 0., 0., 0., 0., 0.
        };
    } else {
        return {
                   talpha,
                   talpha,
                   talpha,
                   0.,
                   0.,
                   0.,
        };
    }
}

void
RheoChainMaterial :: computeDiscreteRelaxationFunction(FloatArray &answer,
                                                       const FloatArray &tSteps,
                                                       double t0, double tr,
                                                       GaussPoint *gp, TimeStep *tStep) const
{
    int size = tSteps.giveSize();
    int nsteps = size;
    FloatArray deltaSigma(size);
    answer.resize(size);
    answer.zero();

    double Jtrt0 = this->computeCreepFunction(t0 + tSteps.at(1), t0, gp, tStep);
    double sig0  = 1. / Jtrt0;
    answer.at(1) = sig0;
    int si = 2;

    double totalDeltaSigma = 0.;
    for ( int k = si; k <= nsteps; k++ ) {
        double sum = 0.;
        for ( int i = si; i <= k - 1; i++ ) {
            double taui;
            if ( i == 1 ) {
                // ??? it seems that tr plays no role because si=2 and the case
                // i=1 or k=1 can never occur
                taui = 0.5 * ( t0 + tSteps.at(i) + tr );
            } else {
                taui = t0 + 0.5 * ( tSteps.at(i) + tSteps.at(i - 1) );
            }

            sum += deltaSigma.at(i) * this->computeCreepFunction(t0 + tSteps.at(k), taui, gp, tStep);
        }

        double tauk;
        if ( k == 1 ) {
            // ??? it seems that tr plays no role because si=2 and the case
            // i=1 or k=1 can never occur
            tauk = 0.5 * ( t0 + tSteps.at(k) + tr );
        } else {
            tauk = t0 + 0.5 * ( tSteps.at(k) + tSteps.at(k - 1) );
        }

        deltaSigma.at(k) = ( sig0 * ( this->computeCreepFunction(t0 + tSteps.at(k), t0, gp, tStep) - Jtrt0 )
                             - sum ) / this->computeCreepFunction(t0 + tSteps.at(k), tauk, gp, tStep);

        totalDeltaSigma += deltaSigma.at(k);
        answer.at(k) = sig0 - totalDeltaSigma;
    }
}


FloatArray
RheoChainMaterial :: generateLogTimeScale(double from, double to, int nsteps)
{
    FloatArray answer(nsteps);
    double help = log(to / from) / nsteps;
    for ( int i = 1; i <= nsteps; i++ ) {
        answer.at(i) = exp(i * help) * from;
    }
    return answer;
}


const FloatArray &
RheoChainMaterial :: giveDiscreteTimes() const
{
    return discreteTimeScale;
}


void
RheoChainMaterial :: giveUnitStiffnessMatrix(FloatMatrix &answer,
                                             GaussPoint *gp,
                                             TimeStep *tStep) const
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
    this->linearElasticMaterial->giveStiffnessMatrix(answer, ElasticStiffness, gp, tStep);
}

void
RheoChainMaterial :: giveUnitComplianceMatrix(FloatMatrix &answer,
                                              GaussPoint *gp,
                                              TimeStep *tStep) const
/*
 * Returns the compliance matrix of an isotropic linear elastic material
 * with the given Poisson ratio (assumed to be unaffected by creep)
 * and a unit value of Young's modulus.
 */
{
    FloatMatrix tangent;
    this->linearElasticMaterial->giveStiffnessMatrix(tangent, ElasticStiffness, gp, tStep);
    answer.beInverseOf(tangent);
}



double
RheoChainMaterial :: giveEparModulus(int iChain) const
{
    /* returns the modulus of unit number iChain,
     * previously computed by updateEparModuli() function
     */
    return EparVal.at(iChain);
}


void
RheoChainMaterial :: updateEparModuli(double tPrime, GaussPoint *gp, TimeStep *tStep) const
{
    /*
     * Computes moduli of individual units in the chain that provide
     * the best approximation of the relaxation or creep function,
     * depending on whether a Maxwell or Kelvin chain is used.
     *
     * INPUTS:
     *
     * tStep - age of material when load is applied ???
     *
     * DESCRIPTION:
     * We store the computed values because they will be used by other material points in subsequent
     * calculations. Their computation is very costly.
     * For a given active time step, there should be requests only for Epar values at time
     * (currentTime+prevTime)*0.5. The previous values are not needed.
     *
     */
    // compute new values and store them in a temporary array for further use
    if ( fabs(tPrime - this->EparValTime) > TIME_DIFF ) {
        this->EparVal = this->computeCharCoefficients(tPrime < 0 ? 1.e-3 : tPrime, gp, tStep);
        this->EparValTime = tPrime;
    }
}


void
RheoChainMaterial :: computeTrueStressIndependentStrainVector(FloatArray &answer,
                                                              GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const
//
// computes the strain due to temperature and shrinkage effects
// (it is called "true" because the "stress-independent strain"
//  will also contain the effect of creep)
//
{
    if ( !Material :: isActivated(tStep) ) {
        answer.zero();
        return;
    }

    // shrinkage strain
    this->giveShrinkageStrainVector(answer, gp, tStep, mode);

    // thermally induced strain
    auto e0 = StructuralMaterial :: computeStressIndependentStrainVector(gp, tStep, mode);
    answer.add(e0);

    if ( e0.giveSize() ) {
#ifdef keep_track_of_strains
        RheoChainMaterialStatus *status = static_cast< RheoChainMaterialStatus * >( this->giveStatus(gp) );
        status->setTempThermalStrain(status->giveThermalStrain()  + e0.at(1) );
#endif
    }
}


FloatArray
RheoChainMaterial :: computeStressIndependentStrainVector(GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const
//
// computes the strain due to temperature, shrinkage and creep effects
// (it is the strain which would occur at the end of the step if the stress
//  during the step was held constant, so it is not really the part of strain
//  independent of stress, but independent of the stress increment)
//
// takes into account the form of the load vector assumed by engngModel (Incremental or Total Load form)
//
{
    // strain due to temperature changes and shrinkage
    FloatArray answer;
    this->computeTrueStressIndependentStrainVector(answer, gp, tStep, mode);
    // strain due to creep
    if ( Material :: isActivated(tStep) ) {
        FloatArray e0;
        this->giveEigenStrainVector(e0, gp, tStep, mode);
        if ( e0.giveSize() ) {
            answer.add(e0);
        }
    }
    return answer;
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
     * log scale. It is necessary to check the condition stepMultiplier <= 10, where Tau(k) = stepMultiplier Tau(k-1)
     *
     *
     */

    int size, nsteps;
    double endTime, Taun1, Tau1, help;

    double stepMultiplier = 10.;

    endTime = this->giveEndOfTimeOfInterest() + relMatAge;
    Taun1 = 0.75 * endTime;

    if ( this->begOfTimeOfInterest == -1 ) {
        this->begOfTimeOfInterest = 0.1;         //default value
    }

    Tau1 = begOfTimeOfInterest;

    if ( Tau1 <= 0 ) {
        OOFEM_ERROR("begOfTimeOfInterest must be a positive number");
    }

    nsteps = ( int ) ( ( log(Taun1) - log(Tau1) ) / log(stepMultiplier) + 1. );
    if ( nsteps < 8 ) {
        nsteps = 8;
    }

    this->nUnits = size = nsteps + 2;

    this->charTimes.resize(size);
    this->charTimes.zero();

    charTimes.at(1) = Tau1;
    charTimes.at(size) = 1.0e10;
    help = ( log(Taun1) - log(Tau1) ) / ( double ) nsteps;
    for ( int i = 1; i <= nsteps; i++ ) {
        charTimes.at(i + 1) = exp(log(Tau1) + help * i);
    }
}

/*
 * void
 * RheoChainMaterial :: giveStiffnessMatrix(FloatMatrix &answer,
 *                                        MatResponseMode rMode,
 *                                        GaussPoint *gp, TimeStep *tStep)
 * {
 *  //
 *  // Returns the incremental material stiffness matrix of the receiver
 *  // PH: the "giveEModulus must be called before any method on elastic material
 *  // otherwise the status would refer to the elastic material and not to the
 *  // viscoelastic one
 *  // in my opinion ElasticStiffness should return incremental stiffness and not unit stiffness
 *  // for this purpose use giveUnitStiffnessMatrix
 *  //
 *  this->giveStatus(gp); // Ensures correct status creation
 *  if ( ( !Material :: isActivated(tStep) ) && ( preCastingTimeMat > 0 ) ) {
 *    auto sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(preCastingTimeMat) );
 *    sMat->giveStiffnessMatrix(answer, rMode, gp, tStep);
 *    // return sMat->give3dMaterialStiffnessMatrix(mode, gp, tStep);
 *  } else {
 *    this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, rMode, gp, tStep);
 *    double incrStiffness = this->giveEModulus(gp, tStep);
 *    answer.times(incrStiffness);
 *    //return incrStiffness * this->linearElasticMaterial->give3dMaterialStiffnessMatrix(mode, gp, tStep);
 *  }
 * }
 */


FloatMatrixF< 6, 6 >
RheoChainMaterial :: give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    //
    // Returns the incremental material stiffness matrix of the receiver
    // PH: the "giveEModulus must be called before any method on elastic material
    // otherwise the status would refer to the elastic material and not to the
    // viscoelastic one
    // in my opinion ElasticStiffness should return incremental stiffness and not unit stiffness
    // for this purpose use giveUnitStiffnessMatrix
    //
    this->giveStatus(gp); // Ensures correct status creation
    if ( ( !Material :: isActivated(tStep) ) && ( preCastingTimeMat > 0 ) ) {
        auto sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(preCastingTimeMat) );
        return sMat->give3dMaterialStiffnessMatrix(mode, gp, tStep);
    } else {
        double incrStiffness = this->giveEModulus(gp, tStep);
        return incrStiffness * this->linearElasticMaterial->give3dMaterialStiffnessMatrix(mode, gp, tStep);
    }
}



FloatMatrixF< 3, 3 >
RheoChainMaterial :: givePlaneStressStiffMtrx(MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *tStep) const
{
    this->giveStatus(gp); // Ensures correct status creation
    if ( ( !Material :: isActivated(tStep) ) && ( preCastingTimeMat > 0 ) ) {
        auto sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(preCastingTimeMat) );
        return sMat->givePlaneStressStiffMtrx(mode, gp, tStep);
    } else {
        double incrStiffness = this->giveEModulus(gp, tStep);
        return incrStiffness * this->linearElasticMaterial->givePlaneStressStiffMtrx(mode, gp, tStep);
    }
}



FloatMatrixF< 4, 4 >
RheoChainMaterial :: givePlaneStrainStiffMtrx(MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *tStep) const
{
    this->giveStatus(gp); // Ensures correct status creation
    if ( ( !Material :: isActivated(tStep) ) && ( preCastingTimeMat > 0 ) ) {
        auto sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(preCastingTimeMat) );
        return sMat->givePlaneStrainStiffMtrx(mode, gp, tStep);
    } else {
        double incrStiffness = this->giveEModulus(gp, tStep);
        return incrStiffness * this->linearElasticMaterial->givePlaneStrainStiffMtrx(mode, gp, tStep);
    }
}



FloatMatrixF< 1, 1 >
RheoChainMaterial :: give1dStressStiffMtrx(MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep) const
{
    this->giveStatus(gp); // Ensures correct status creation
    if ( ( !Material :: isActivated(tStep) ) && ( preCastingTimeMat > 0 ) ) {
        auto sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(preCastingTimeMat) );
        return sMat->give1dStressStiffMtrx(mode, gp, tStep);
    } else {
        double incrStiffness = this->giveEModulus(gp, tStep);
        return incrStiffness * this->linearElasticMaterial->give1dStressStiffMtrx(mode, gp, tStep);
    }
}

MaterialStatus *
RheoChainMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new RheoChainMaterialStatus(gp, nUnits);
}


void
RheoChainMaterial :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);

    // if the casting time is not defined, we set it to zero such that it concides
    // with the beginning of the computational time
    if ( !ir.hasField(_IFT_Material_castingtime) ) {
        this->castingTime = 0.;
    }

    this->talpha = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, talpha, _IFT_RheoChainMaterial_talpha);

    if ( ir.hasField(_IFT_RheoChainMaterial_lattice) ) {
        lattice = true;
        this->alphaOne = 1.;
        this->alphaTwo = 1.;
        IR_GIVE_OPTIONAL_FIELD(ir, alphaOne, _IFT_RheoChainMaterial_alphaOne);
        IR_GIVE_OPTIONAL_FIELD(ir, alphaTwo, _IFT_RheoChainMaterial_alphaTwo);
    } else {
        lattice = false;
        IR_GIVE_FIELD(ir, nu, _IFT_RheoChainMaterial_n);
    }

    IR_GIVE_FIELD(ir, relMatAge, _IFT_RheoChainMaterial_relmatage);
    this->begOfTimeOfInterest = -1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, begOfTimeOfInterest, _IFT_RheoChainMaterial_begoftimeofinterest);
    this->endOfTimeOfInterest = -1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, endOfTimeOfInterest, _IFT_RheoChainMaterial_endoftimeofinterest);
    IR_GIVE_FIELD(ir, timeFactor, _IFT_RheoChainMaterial_timefactor); // solution time/timeFactor should give time in days

    // sets up nUnits variable and characteristic times array (retardation/relaxation times)
    this->computeCharTimes();

    // sets up discrete times
    double endTime = this->giveEndOfTimeOfInterest();
    this->discreteTimeScale = this->generateLogTimeScale(this->begOfTimeOfInterest, endTime, MNC_NPOINTS - 1);

    this->giveLinearElasticMaterial();
    ///@warning Stiffness is time dependant, so the variable changes with time.
    //ph !!! why was it put here?
    //this->updateEparModuli(0.); // stiffnesses are time independent (evaluated at time t = 0.)
}

//LinearElasticMaterial *
StructuralMaterial *
RheoChainMaterial :: giveLinearElasticMaterial()
{
    if ( linearElasticMaterial == NULL ) {
        if ( this->lattice ) {
            linearElasticMaterial = new LatticeLinearElastic(this->giveNumber(),
                                                             this->giveDomain(),
                                                             1.0, this->alphaOne, this->alphaTwo);
        } else {
            linearElasticMaterial = new IsotropicLinearElasticMaterial(this->giveNumber(),
                                                                       this->giveDomain(),
                                                                       1.0, this->nu);
        }
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

void
RheoChainMaterial :: saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp)
{
    Material :: saveIPContext(stream, mode, gp);
    giveLinearElasticMaterial()->saveIPContext(stream, mode, gp);
}

void
RheoChainMaterial :: restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp)
{
    Material :: restoreIPContext(stream, mode, gp);
    // invoke possible restoring of statuses for submaterials
    giveLinearElasticMaterial()->restoreIPContext(stream, mode, gp);
}


int
RheoChainMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_ThermalStrainTensor ) {
        RheoChainMaterialStatus *status = static_cast< RheoChainMaterialStatus * >( this->giveStatus(gp) );
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveThermalStrain();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }

    // return 1; // to make the compiler happy
}

/****************************************************************************************/

RheoChainMaterialStatus :: RheoChainMaterialStatus(GaussPoint *g, int nunits) :
    StructuralMaterialStatus(g),
    nUnits(nunits),
    hiddenVars(nUnits),
    tempHiddenVars(nUnits)
{}


void
RheoChainMaterialStatus :: letTempHiddenVarsVectorBe(int i, FloatArray &valueArray)
{
    // Sets the i:th hidden variables vector to valueArray.
#ifdef DEBUG
    if ( i > nUnits ) {
        OOFEM_ERROR("unit number exceeds the specified limit");
    }
#endif
    this->tempHiddenVars [ i - 1 ] = valueArray;
}

void
RheoChainMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    // printing of hidden variables

    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "{hidden variables: ");
    for ( int i = 0; i < nUnits; i++ ) {
        fprintf(file, "{ ");
        for ( auto &val : hiddenVars [ i ] ) {
            fprintf(file, "%f ", val);
        }

        fprintf(file, "} ");
    }

    if ( shrinkageStrain.isNotEmpty() ) {
        fprintf(file, "shrinkageStrain: {");
        for ( auto &val : shrinkageStrain ) {
            fprintf(file, "%f ", val);
        }

        fprintf(file, "} ");
    }

    fprintf(file, "}\n");
}


void
RheoChainMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);

    for ( int i = 0; i < nUnits; i++ ) {
        this->hiddenVars [ i ] = this->tempHiddenVars [ i ];
    }

    currentTime = -1.e20;

#ifdef keep_track_of_strains
    thermalStrain = tempThermalStrain;
    tempThermalStrain = 0.;
#endif
}

void
RheoChainMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();

#ifdef keep_track_of_strains
    tempThermalStrain = 0.;
#endif
}

void
RheoChainMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;
    for ( int i = 0; i < nUnits; i++ ) {
        if ( ( iores = hiddenVars [ i ].storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( ( iores = shrinkageStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


void
RheoChainMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus :: restoreContext(stream, mode);

    contextIOResultType iores;
    for ( int i = 0; i < nUnits; i++ ) {
        if ( ( iores = hiddenVars [ i ].restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( ( iores = shrinkageStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}
} // end namespace oofem
