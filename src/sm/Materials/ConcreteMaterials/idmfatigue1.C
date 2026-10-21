/*
 * Fatigue extension of IsotropicDamageMaterial1. See idmfatigue1.h for the model description.
 * Derived from the OOFEM framework (LGPL 2.1+), same licence as idm1.
 */

#include "idmfatigue1.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "domain.h"
#include "function.h"
#include "timestep.h"

#include <algorithm>
#include <cmath>

namespace oofem {
REGISTER_Material(IsotropicDamageMaterialFatigue1);

/////////////////////////////////////////////////////////////////////////////////////////////////
// Fatigue law (shared with the nonlocal material)
/////////////////////////////////////////////////////////////////////////////////////////////////

void
IDMFatigueLaw :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    IR_GIVE_FIELD(ir, A, _IFT_IsotropicDamageMaterialFatigue1_A);
    IR_GIVE_OPTIONAL_FIELD(ir, m, _IFT_IsotropicDamageMaterialFatigue1_m);
    IR_GIVE_OPTIONAL_FIELD(ir, s0, _IFT_IsotropicDamageMaterialFatigue1_s0);
    IR_GIVE_OPTIONAL_FIELD(ir, kd, _IFT_IsotropicDamageMaterialFatigue1_kd);
    IR_GIVE_OPTIONAL_FIELD(ir, dmax, _IFT_IsotropicDamageMaterialFatigue1_dmax);
    IR_GIVE_OPTIONAL_FIELD(ir, cj, _IFT_IsotropicDamageMaterialFatigue1_cj);
    IR_GIVE_OPTIONAL_FIELD(ir, cjFunc, _IFT_IsotropicDamageMaterialFatigue1_cjfunc);

    if ( A < 0. ) {
        throw ValueInputException(ir, _IFT_IsotropicDamageMaterialFatigue1_A, "must be non-negative");
    }
    if ( m < 0. ) {
        throw ValueInputException(ir, _IFT_IsotropicDamageMaterialFatigue1_m, "must be non-negative");
    }
    if ( s0 < 0. || s0 >= 1. ) {
        throw ValueInputException(ir, _IFT_IsotropicDamageMaterialFatigue1_s0, "must satisfy 0 <= s0 < 1");
    }
    if ( kd < 0. ) {
        throw ValueInputException(ir, _IFT_IsotropicDamageMaterialFatigue1_kd, "must be non-negative");
    }
    if ( dmax <= 0. || dmax >= 1. ) {
        throw ValueInputException(ir, _IFT_IsotropicDamageMaterialFatigue1_dmax, "must satisfy 0 < fdmax < 1");
    }
    if ( cj < 1. ) {
        throw ValueInputException(ir, _IFT_IsotropicDamageMaterialFatigue1_cj, "must be >= 1");
    }
}

void
IDMFatigueLaw :: giveInputRecord(DynamicInputRecord &input) const
{
    input.setField(A, _IFT_IsotropicDamageMaterialFatigue1_A);
    input.setField(m, _IFT_IsotropicDamageMaterialFatigue1_m);
    input.setField(s0, _IFT_IsotropicDamageMaterialFatigue1_s0);
    input.setField(kd, _IFT_IsotropicDamageMaterialFatigue1_kd);
    input.setField(dmax, _IFT_IsotropicDamageMaterialFatigue1_dmax);
    input.setField(cj, _IFT_IsotropicDamageMaterialFatigue1_cj);
    input.setField(cjFunc, _IFT_IsotropicDamageMaterialFatigue1_cjfunc);
}

double
IDMFatigueLaw :: primitive(double x) const
{
    return std::pow(std::max(x - s0, 0.), m + 1.);
}

double
IDMFatigueLaw :: currentCj(Domain *d, TimeStep *tStep) const
{
    double f = cj;
    if ( cjFunc > 0 ) {
        f *= d->giveFunction(cjFunc)->evaluateAtTime( tStep->giveIntrinsicTime() );
    }
    return std::max(f, 1.);
}

double
IDMFatigueLaw :: advance(double D_old, double lo, double hi, double e0, double cjNow) const
{
    // Rate law   dD * (1-D)^kd = A * dF(x),  F(x) = <x - s0>^(m+1),  x = eq / e0
    // is separable, so it is integrated exactly:
    // (1-D_new)^(kd+1) = (1-D_old)^(kd+1) - (kd+1) * A * cjNow * [F(x_hi) - F(x_lo)]
    if ( !( hi > lo ) || A <= 0. ) {
        return D_old;
    }
    const double dPhi = A * cjNow * ( primitive(hi / e0) - primitive(lo / e0) );
    if ( dPhi <= 0. ) {
        return D_old;
    }
    const double p = kd + 1.;
    const double arg = std::pow(1. - D_old, p) - p * dPhi;
    const double argMin = std::pow(1. - dmax, p);
    const double D = 1. - std::pow(std::max(arg, argMin), 1. / p);
    return std::max(D, D_old); // damage never decreases
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Status
/////////////////////////////////////////////////////////////////////////////////////////////////

void
IsotropicDamageMaterialFatigue1Status :: initTempStatus()
{
    IsotropicDamageMaterial1Status :: initTempStatus();
    tempEqStrainPrev = eqStrainPrev;
    tempStaticDamage = staticDamage;
}

void
IsotropicDamageMaterialFatigue1Status :: updateYourself(TimeStep *tStep)
{
    IsotropicDamageMaterial1Status :: updateYourself(tStep);
    eqStrainPrev = tempEqStrainPrev;
    staticDamage = tempStaticDamage;
}

void
IsotropicDamageMaterialFatigue1Status :: saveContext(DataStream &stream, ContextMode mode)
{
    IsotropicDamageMaterial1Status :: saveContext(stream, mode);
    if ( !stream.write(& eqStrainPrev, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    if ( !stream.write(& staticDamage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}

void
IsotropicDamageMaterialFatigue1Status :: restoreContext(DataStream &stream, ContextMode mode)
{
    IsotropicDamageMaterial1Status :: restoreContext(stream, mode);
    if ( !stream.read(& eqStrainPrev, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    if ( !stream.read(& staticDamage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    tempEqStrainPrev = eqStrainPrev;
    tempStaticDamage = staticDamage;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Material
/////////////////////////////////////////////////////////////////////////////////////////////////

IsotropicDamageMaterialFatigue1 :: IsotropicDamageMaterialFatigue1(int n, Domain *d) :
    IsotropicDamageMaterial1(n, d)
{ }

void
IsotropicDamageMaterialFatigue1 :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    IsotropicDamageMaterial1 :: initializeFrom(ir);
    law.initializeFrom(ir);
}

void
IsotropicDamageMaterialFatigue1 :: giveInputRecord(DynamicInputRecord &input)
{
    IsotropicDamageMaterial1 :: giveInputRecord(input);
    law.giveInputRecord(input);
}

std::unique_ptr<MaterialStatus>
IsotropicDamageMaterialFatigue1 :: CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<IsotropicDamageMaterialFatigue1Status>(gp);
}

void
IsotropicDamageMaterialFatigue1 :: giveFatigueStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep)
{
    auto status = static_cast< IsotropicDamageMaterialFatigue1Status * >( this->giveStatus(gp) );
    status->initTempStatus();

    // mechanical (stress-producing) part of the strain, i.e. without thermal/shrinkage parts
    FloatArray strainVector;
    this->giveStressDependentPartOfStrainVector(strainVector, gp, totalStrain, tStep, VM_Total);

    const double e0 = this->give(e0_ID, gp);
    const double eq = this->computeEquivalentStrain(strainVector, gp, tStep);

    // converged state at the beginning of the step (never modified inside the equilibrium iterations)
    const double D_old = status->giveDamage();
    const double kap_old = status->giveKappa();
    const double eq_old = status->giveEqStrainPrev();
    // damage envelope; kappa itself starts from 0, but nothing can be "static" below e0
    const double env = std::max(kap_old, e0);

    // (1) FATIGUE: part of the path that is increasing and lies below the envelope
    double D = law.advance(D_old, eq_old, std::min(eq, env), e0, law.currentCj(this->giveDomain(), tStep));

    // (2) STATIC: the equivalent strain sets a new maximum -> usual idm1 damage law, added as an increment
    const double kap = std::max(kap_old, eq);
    double Ds = status->giveStaticDamage();
    // crack direction and element size are stored in the status by initDamaged (idm1 does this at
    // damage onset, or when Le is not set yet); also needed when damage is caused by fatigue only
    if ( D > D_old || kap > env ) {
        this->initDamaged(std::max(kap, e0 * ( 1. + 1.e-9 )), strainVector, gp);
    }
    if ( kap > env ) {
        const double g = this->computeDamageParam(kap, strainVector, gp);
        if ( g > Ds ) {
            D += g - Ds;
            Ds = g;
        }
    }
    D = std::min(D, law.dmax);

    // (3) stress, secant stiffness
    FloatMatrix de;
    this->linearElasticMaterial->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
    de.times(1. - D);
    answer.beProductOf(de, strainVector);

    // (4) store temporary state
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->setTempKappa(kap);
    status->setTempDamage(D);
    status->setTempEqStrainPrev(eq);
    status->setTempStaticDamage(Ds);
}

// Mode-specific entry points: same forwarding pattern as in IsotropicDamageMaterial.
FloatArrayF<6>
IsotropicDamageMaterialFatigue1 :: giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray answer;
    const_cast< IsotropicDamageMaterialFatigue1 * >( this )->giveFatigueStressVector(answer, gp, strain, tStep);
    return FloatArrayF<6>(answer);
}

FloatArrayF<4>
IsotropicDamageMaterialFatigue1 :: giveRealStressVector_PlaneStrain(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray answer;
    const_cast< IsotropicDamageMaterialFatigue1 * >( this )->giveFatigueStressVector(answer, gp, strain, tStep);
    return FloatArrayF<4>(answer);
}

FloatArray
IsotropicDamageMaterialFatigue1 :: giveRealStressVector_StressControl(const FloatArray &strain, const IntArray &strainControl, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray answer;
    const_cast< IsotropicDamageMaterialFatigue1 * >( this )->giveFatigueStressVector(answer, gp, strain, tStep);
    return answer;
}

FloatArrayF<3>
IsotropicDamageMaterialFatigue1 :: giveRealStressVector_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray answer;
    const_cast< IsotropicDamageMaterialFatigue1 * >( this )->giveFatigueStressVector(answer, gp, strain, tStep);
    return FloatArrayF<3>(answer);
}

FloatArrayF<1>
IsotropicDamageMaterialFatigue1 :: giveRealStressVector_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray answer;
    const_cast< IsotropicDamageMaterialFatigue1 * >( this )->giveFatigueStressVector(answer, gp, strain, tStep);
    return FloatArrayF<1>(answer);
}
} // end namespace oofem
