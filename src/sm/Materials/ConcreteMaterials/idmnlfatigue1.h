/*
 * Nonlocal isotropic damage model with fatigue (IDNLMaterial + fatigue law of idmfatigue1.h).
 *
 * The fatigue and static damage are driven by the NONLOCAL equivalent strain (averaged as in idmnl1),
 * so both the static softening and the fatigue growth are regularised by the same internal length.
 * Only the standard formulation (averaging of equivalent strain) is supported; averaging of
 * compliance or damage is rejected at input, because the fatigue law needs a strain measure.
 *
 * Input: all fields of idmnl1 (including the nonlocal ones, e.g. r) plus fa, fm, fs0, fkd, fdmax, fcj
 * (see idmfatigue1.h). Record name: idmnlfat1.
 */

#ifndef idmnlfatigue1_h
#define idmnlfatigue1_h

#include "idmnl1.h"
#include "idmfatigue1.h"

///@name Input fields for IDNLFatigueMaterial1
//@{
#define _IFT_IDNLFatigueMaterial1_Name "idmnlfat1"
//@}

namespace oofem {
/// Status of IDNLFatigueMaterial1: IDNL status + fatigue history variables.
class IDNLFatigueMaterial1Status : public IDNLMaterialStatus
{
protected:
    double eqStrainPrev = 0., tempEqStrainPrev = 0.;
    double staticDamage = 0., tempStaticDamage = 0.;

public:
    IDNLFatigueMaterial1Status(GaussPoint *g) : IDNLMaterialStatus(g) { }

    const char *giveClassName() const override { return "IDNLFatigueMaterial1Status"; }

    double giveEqStrainPrev() const { return eqStrainPrev; }
    void setTempEqStrainPrev(double v) { tempEqStrainPrev = v; }
    double giveStaticDamage() const { return staticDamage; }
    void setTempStaticDamage(double v) { tempStaticDamage = v; }
    double giveFatigueDamage() const { double d = this->giveDamage() - staticDamage; return d > 0. ? d : 0.; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};

/// Nonlocal isotropic damage model with fatigue.
class IDNLFatigueMaterial1 : public IDNLMaterial
{
protected:
    IDMFatigueLaw law;

public:
    IDNLFatigueMaterial1(int n, Domain *d);

    const char *giveClassName() const override { return "IDNLFatigueMaterial1"; }
    const char *giveInputRecordName() const override { return _IFT_IDNLFatigueMaterial1_Name; }
    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    void giveFatigueStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep);
    /// Stores crack direction and angle in the status (IDNLMaterial::initDamaged is empty, so idmnl1 never sets them).
    void setCrackDirection(const FloatArray &strainVector, GaussPoint *gp);

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<4> giveRealStressVector_PlaneStrain(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArray giveRealStressVector_StressControl(const FloatArray &strain, const IntArray &strainControl, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<IDNLFatigueMaterial1Status>(gp); }
};
} // end namespace oofem
#endif // idmnlfatigue1_h
