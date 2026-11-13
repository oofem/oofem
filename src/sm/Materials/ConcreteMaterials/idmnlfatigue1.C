/*
 * Nonlocal isotropic damage model with fatigue. See idmnlfatigue1.h.
 * Derived from the OOFEM framework (LGPL 2.1+), same licence as idm1/idmnl1.
 */

#include "idmnlfatigue1.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "datastream.h"
#include "contextioerr.h"
#include "error.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

#include <algorithm>
#include <cmath>

namespace oofem {
REGISTER_Material(IDNLFatigueMaterial1);

/////////////////////////////////////////////////////////////////////////////////////////////////
// Status
/////////////////////////////////////////////////////////////////////////////////////////////////

void
IDNLFatigueMaterial1Status :: initTempStatus()
{
    IDNLMaterialStatus :: initTempStatus();
    tempEqStrainPrev = eqStrainPrev;
    tempStaticDamage = staticDamage;
}

void
IDNLFatigueMaterial1Status :: updateYourself(TimeStep *tStep)
{
    IDNLMaterialStatus :: updateYourself(tStep);
    eqStrainPrev = tempEqStrainPrev;
    staticDamage = tempStaticDamage;
}

void
IDNLFatigueMaterial1Status :: saveContext(DataStream &stream, ContextMode mode)
{
    IDNLMaterialStatus :: saveContext(stream, mode);
    if ( !stream.write(& eqStrainPrev, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    if ( !stream.write(& staticDamage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}

void
IDNLFatigueMaterial1Status :: restoreContext(DataStream &stream, ContextMode mode)
{
    IDNLMaterialStatus :: restoreContext(stream, mode);
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

IDNLFatigueMaterial1 :: IDNLFatigueMaterial1(int n, Domain *d) :
    IDNLMaterial(n, d)
{ }

void
IDNLFatigueMaterial1 :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    IDNLMaterial :: initializeFrom(ir);
    law.initializeFrom(ir);

    // the fatigue law needs a strain measure -> only averaging of the equivalent strain is allowed
    if ( averagedVar == AVT_Compliance || averagedVar == AVT_Damage ) {
        OOFEM_ERROR("idmnlfat1 supports only averaging of the equivalent strain (not compliance or damage)");
    }
}

void
IDNLFatigueMaterial1 :: giveInputRecord(DynamicInputRecord &input)
{
    IDNLMaterial :: giveInputRecord(input);
    law.giveInputRecord(input);
}

void
IDNLFatigueMaterial1 :: setCrackDirection(const FloatArray &strainVector, GaussPoint *gp)
{
    // same definition as in IsotropicDamageMaterial1::initDamaged (without the Griffith special case)
    auto status = static_cast< IDNLFatigueMaterial1Status * >( this->giveStatus(gp) );
    FloatArray fullStrain, principalStrains, crackVect, crackPlaneNormal;
    FloatMatrix principalDir;
    StructuralMaterial :: giveFullSymVectorForm( fullStrain, strainVector, gp->giveMaterialMode() );
    this->computePrincipalValDir(principalStrains, principalDir, fullStrain, principal_strain);

    // index of the maximum principal strain -> normal to the crack plane
    int indx = 1;
    for ( int i = 2; i <= 3; i++ ) {
        if ( principalStrains.at(i) > principalStrains.at(indx) ) {
            indx = i;
        }
    }
    crackPlaneNormal.beColumnOf(principalDir, indx);

    // minimal non-zero principal strain -> crack direction
    indx = 1;
    for ( int i = 2; i <= 3; i++ ) {
        if ( principalStrains.at(i) < principalStrains.at(indx) && std::fabs( principalStrains.at(i) ) > 1.e-10 ) {
            indx = i;
        }
    }
    crackVect.beColumnOf(principalDir, indx);
    status->setCrackVector(crackVect);

    double ca = M_PI / 2.;
    if ( crackPlaneNormal.at(1) != 0.0 ) {
        ca = std::atan( crackPlaneNormal.at(2) / crackPlaneNormal.at(1) );
    }
    status->setCrackAngle(ca);
}

void
IDNLFatigueMaterial1 :: giveFatigueStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep)
{
    auto status = static_cast< IDNLFatigueMaterial1Status * >( this->giveStatus(gp) );
    status->initTempStatus();

    // mechanical (stress-producing) part of the strain, i.e. without thermal/shrinkage parts
    FloatArray strainVector;
    this->giveStressDependentPartOfStrainVector(strainVector, gp, totalStrain, tStep, VM_Total);

    const double e0 = this->give(e0_ID, gp);
    // NONLOCAL equivalent strain (IDNLMaterial::computeEquivalentStrain averages over the neighbourhood)
    const double eq = this->computeEquivalentStrain(strainVector, gp, tStep);

    const double D_old = status->giveDamage();
    const double kap_old = status->giveKappa();
    const double eq_old = status->giveEqStrainPrev();
    const double env = std::max(kap_old, e0);

    // (1) fatigue: increasing part of the nonlocal path below the envelope
    double D = law.advance(D_old, eq_old, std::min(eq, env), e0, law.currentCj(this->FEMComponent :: giveDomain(), tStep));

    // (2) static: new maximum of the nonlocal equivalent strain (initDamaged is a no-op in IDNLMaterial)
    const double kap = std::max(kap_old, eq);
    double Ds = status->giveStaticDamage();
    if ( kap > env ) {
        const double g = this->computeDamageParam(kap, strainVector, gp);
        if ( g > Ds ) {
            D += g - Ds;
            Ds = g;
        }
    }
    D = std::min(D, law.dmax);

    // crack direction at damage onset (static or fatigue)
    if ( D_old == 0. && D > 0. ) {
        this->setCrackDirection(strainVector, gp);
    }

    // (3) stress, secant stiffness
    FloatMatrix de;
    this->linearElasticMaterial->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
    de.times(1. - D);
    answer.beProductOf(de, strainVector);

    // (4) temporary state
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->setTempKappa(kap);
    status->setTempDamage(D);
    status->setTempEqStrainPrev(eq);
    status->setTempStaticDamage(Ds);
}

FloatArrayF<6>
IDNLFatigueMaterial1 :: giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray answer;
    const_cast< IDNLFatigueMaterial1 * >( this )->giveFatigueStressVector(answer, gp, strain, tStep);
    return FloatArrayF<6>(answer);
}

FloatArrayF<4>
IDNLFatigueMaterial1 :: giveRealStressVector_PlaneStrain(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray answer;
    const_cast< IDNLFatigueMaterial1 * >( this )->giveFatigueStressVector(answer, gp, strain, tStep);
    return FloatArrayF<4>(answer);
}

FloatArray
IDNLFatigueMaterial1 :: giveRealStressVector_StressControl(const FloatArray &strain, const IntArray &strainControl, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray answer;
    const_cast< IDNLFatigueMaterial1 * >( this )->giveFatigueStressVector(answer, gp, strain, tStep);
    return answer;
}

FloatArrayF<3>
IDNLFatigueMaterial1 :: giveRealStressVector_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray answer;
    const_cast< IDNLFatigueMaterial1 * >( this )->giveFatigueStressVector(answer, gp, strain, tStep);
    return FloatArrayF<3>(answer);
}

FloatArrayF<1>
IDNLFatigueMaterial1 :: giveRealStressVector_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray answer;
    const_cast< IDNLFatigueMaterial1 * >( this )->giveFatigueStressVector(answer, gp, strain, tStep);
    return FloatArrayF<1>(answer);
}
} // end namespace oofem
