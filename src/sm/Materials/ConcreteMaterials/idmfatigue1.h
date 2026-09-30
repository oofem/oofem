/*
 * Fatigue extension of the OOFEM isotropic damage model IsotropicDamageMaterial1 (idm1).
 *
 * Model (threshold / loading-surface form with a fatigue term):
 *   - sigma = (1 - D) * De : eps                       (secant, scalar isotropic damage)
 *   - kappa = max equivalent strain so far (as in idm1), envelope = max(kappa, e0)
 *   - STATIC damage grows when the equivalent strain sets a new maximum (any idm1 damage law)
 *   - FATIGUE damage grows when the equivalent strain is INCREASING but still BELOW the envelope:
 *
 *        dD_f * (1 - D)^kd = A * d( <eq/e0 - s0>^(m+1) )
 *
 *     integrated exactly over each step (see IDMFatigueLaw::advance), so the result is independent
 *     of step size.
 *
 * Additional input fields (shared by the local "idmfat1" and the nonlocal "idmnlfat1" material):
 *   fa     (required) fatigue amplitude A [-]
 *   fm     exponent m (default 2)
 *   fs0    fatigue threshold s0 (default 0.5), must satisfy 0 <= s0 < 1
 *   fkd    acceleration exponent kd (default 1)
 *   fdmax  cap on total damage (default 0.999)
 *   fcj    cycle-scaling factor (default 1). Each simulated cycle counts as fcj cycles. APPROXIMATE
 *          cycle jump, valid only while stress redistribution between jumps is negligible.
 *   fcjfunc  optional number of a time function f(t); the factor used is fcj*f(t) (at least 1),
 *          so the jump can be changed during the analysis. Real cycles = integral of the factor.
 */

#ifndef idmfatigue1_h
#define idmfatigue1_h

#include "idm1.h"

///@name Input fields for the fatigue law and for IsotropicDamageMaterialFatigue1
//@{
#define _IFT_IsotropicDamageMaterialFatigue1_Name "idmfat1"
#define _IFT_IsotropicDamageMaterialFatigue1_A "fa"
#define _IFT_IsotropicDamageMaterialFatigue1_m "fm"
#define _IFT_IsotropicDamageMaterialFatigue1_s0 "fs0"
#define _IFT_IsotropicDamageMaterialFatigue1_kd "fkd"
#define _IFT_IsotropicDamageMaterialFatigue1_dmax "fdmax"
#define _IFT_IsotropicDamageMaterialFatigue1_cj "fcj"
#define _IFT_IsotropicDamageMaterialFatigue1_cjfunc "fcjfunc"
//@}

namespace oofem {
/**
 * Fatigue damage law shared by the local and nonlocal materials: parameters, input handling and
 * the exact integration of the fatigue increment.
 */
struct IDMFatigueLaw
{
    double A = 0.;        ///< amplitude
    double m = 2.;        ///< exponent
    double s0 = 0.5;      ///< threshold (normalised by e0)
    double kd = 1.;       ///< acceleration exponent
    double dmax = 0.999;  ///< cap on total damage
    double cj = 1.;       ///< cycle scaling factor
    int cjFunc = 0;       ///< optional time function multiplying cj (0 = none)

    void initializeFrom(const std::shared_ptr<InputRecord> &ir);
    void giveInputRecord(DynamicInputRecord &input) const;
    /// Primitive of the fatigue integrand, x = equivalent strain / e0.
    double primitive(double x) const;
    /**
     * Damage after a fatigue increment along the increasing strain path lo -> hi (both below the
     * envelope). Exact solution of dD (1-D)^kd = A dF(x); returns D_old if nothing accumulates.
     */
    double advance(double D_old, double lo, double hi, double e0, double cjNow) const;
    /// Current cycle-scaling factor: cj * f(t), where f is the time function fcjfunc (if given); at least 1.
    double currentCj(Domain *d, TimeStep *tStep) const;
};

/**
 * Status of IsotropicDamageMaterialFatigue1. Adds two history variables to idm1 status:
 *  - equivalent strain at the end of the previous converged step (to know the loading path),
 *  - static damage g(kappa) (envelope value), so that D - g(kappa) is the fatigue part.
 */
class IsotropicDamageMaterialFatigue1Status : public IsotropicDamageMaterial1Status
{
protected:
    double eqStrainPrev = 0., tempEqStrainPrev = 0.;
    double staticDamage = 0., tempStaticDamage = 0.;

public:
    IsotropicDamageMaterialFatigue1Status(GaussPoint *g) : IsotropicDamageMaterial1Status(g) { }

    const char *giveClassName() const override { return "IsotropicDamageMaterialFatigue1Status"; }

    double giveEqStrainPrev() const { return eqStrainPrev; }
    void setTempEqStrainPrev(double v) { tempEqStrainPrev = v; }
    double giveStaticDamage() const { return staticDamage; }
    void setTempStaticDamage(double v) { tempStaticDamage = v; }
    /// Fatigue part of the damage (converged): D - g(kappa), never negative.
    double giveFatigueDamage() const { double d = this->giveDamage() - staticDamage; return d > 0. ? d : 0.; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};

/**
 * Local isotropic damage model for concrete with fatigue (see file header).
 */
class IsotropicDamageMaterialFatigue1 : public IsotropicDamageMaterial1
{
protected:
    IDMFatigueLaw law;

public:
    IsotropicDamageMaterialFatigue1(int n, Domain *d);

    const char *giveClassName() const override { return "IsotropicDamageMaterialFatigue1"; }
    const char *giveInputRecordName() const override { return _IFT_IsotropicDamageMaterialFatigue1_Name; }
    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    /// Common core of all stress evaluations (non-const, like the forwarding target in IsotropicDamageMaterial).
    void giveFatigueStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep);

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<4> giveRealStressVector_PlaneStrain(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArray giveRealStressVector_StressControl(const FloatArray &strain, const IntArray &strainControl, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;
};
} // end namespace oofem
#endif // idmfatigue1_h
