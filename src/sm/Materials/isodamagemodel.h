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

#ifndef isodamagemodel_h
#define isodamagemodel_h

// this turns on or off a bunch of internal variables
// that allow tracing the distribution of dissipated energy
// (can be turned off if such information is not needed)
#define keep_track_of_dissipated_energy

#include "material.h"
#include "sm/Materials/linearelasticmaterial.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"

///@name Input fields for IsotropicDamageMaterial
//@{
#define _IFT_IsotropicDamageMaterial_talpha "talpha"
#define _IFT_IsotropicDamageMaterial_maxOmega "maxomega"
#define _IFT_IsotropicDamageMaterial_permstrain "ps"
//@}

namespace oofem {
class GaussPoint;

/**
 * This class implements associated Material Status to IsotropicDamageMaterial.
 * Stores a scalar damage and hardening variable (and possible extra information).
 */
class IsotropicDamageMaterialStatus : public StructuralMaterialStatus
{
protected:
    /// Scalar measure of the largest strain level ever reached in material.
    double kappa = 0.;
    /// Non-equilibrated scalar measure of the largest strain level.
    double tempKappa = 0.;
    /// Damage level of material.
    double damage = 0.;
    /// Non-equilibrated damage level of material.
    double tempDamage = 0.;
    /**
     * Characteristic element length,
     * computed when damage initialized from direction of
     * maximum positive principal strain. Fixed during further loading.
     */
    double le = 0.;
    /// Angle characterizing the crack direction.
    double crack_angle = -1000.0;
    /// Crack orientation normalized to damage magnitude. This is useful for plotting cracks as a vector field (paraview etc.).
    FloatArrayF<3> crackVector;

#ifdef keep_track_of_dissipated_energy
    /// Density of total work done by stresses on strain increments.
    double stressWork = 0.;
    /// Non-equilibrated density of total work done by stresses on strain increments.
    double tempStressWork = 0.;
    /// Density of dissipated work.
    double dissWork = 0.;
    /// Non-equilibrated density of dissipated work.
    double tempDissWork = 0.;
#endif

public:
    /// Constructor
    IsotropicDamageMaterialStatus(GaussPoint *g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    /// Returns the last equilibrated scalar measure of the largest strain level.
    double giveKappa() const { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level.
    double giveTempKappa() const { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value.
    void setTempKappa(double newKappa) { tempKappa = newKappa; }
    /// Returns the last equilibrated damage level.
    double giveDamage() const { return damage; }
    /// Returns the temp. damage level.
    double giveTempDamage() const { return tempDamage; }
    /// Sets the temp damage level to given value.
    void setTempDamage(double newDamage) { tempDamage = newDamage; }

    /// Returns characteristic length stored in receiver.
    double giveLe() const { return le; }
    /// Sets characteristic length to given value.
    void setLe(double ls) { le = ls; }
    /// Returns crack angle stored in receiver.
    double giveCrackAngle() const { return crack_angle; }
    /// Sets crack angle to given value.
    void setCrackAngle(double ca) { crack_angle = ca; }
    /// Returns crack vector stored in receiver. This is useful for plotting cracks as a vector field (paraview etc.).
    FloatArrayF<3> giveCrackVector() const { return crackVector * damage; }
    /// Sets crack vector to given value. This is useful for plotting cracks as a vector field (paraview etc.).
    void setCrackVector(const FloatArrayF<3> &cv) { crackVector = cv; }

#ifdef keep_track_of_dissipated_energy
    /// Returns the density of total work of stress on strain increments.
    double giveStressWork() const { return stressWork; }
    /// Returns the temp density of total work of stress on strain increments.
    double giveTempStressWork() const { return tempStressWork; }
    /// Sets the density of total work of stress on strain increments to given value.
    void setTempStressWork(double w) { tempStressWork = w; }
    /// Returns the density of dissipated work.
    double giveDissWork() const { return dissWork; }
    /// Returns the density of temp dissipated work.
    double giveTempDissWork() const { return tempDissWork; }
    /// Sets the density of dissipated work to given value.
    void setTempDissWork(double w) { tempDissWork = w; }
    /// Computes the increment of total stress work and of dissipated work.
    void computeWork(GaussPoint *gp);
#endif

    const char *giveClassName() const override { return "IsotropicDamageMaterialModelStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * Base class representing general isotropic damage model.
 * It is based on isotropic damage concept, assuming that damage evolution law
 * is postulated in explicit form, relation damage parameter (omega) to scalar measure
 * of the largest strain level ever reached in material (kappa).
 */
class IsotropicDamageMaterial : public StructuralMaterial
{
protected:
    /// Coefficient of thermal dilatation.
    double tempDillatCoeff = 0.;

    /// Maximum limit on omega. The purpose is elimination of a too compliant material which may cause convergence problems. Set to something like 0.99 if needed.
    double maxOmega = 0.999999;

    /// Indicator of the type of permanent strain formulation (0 = standard damage with no permanent strain)
    int permStrain = 0;

    /// Reference to bulk (undamaged) material
    LinearElasticMaterial *linearElasticMaterial = nullptr;
    /**
     * Variable controlling type of loading/unloading law, default set to idm_strainLevel
     * defines the two two possibilities:
     * - idm_strainLevelCR the unloading takes place, when strain level is smaller than the largest level ever reached;
     * - idm_damageLevelCR the unloading takes place, when damage level is smaller than the largest damage ever  reached;
     */
    enum loaUnloCriterium { idm_strainLevelCR, idm_damageLevelCR } llcriteria = idm_strainLevelCR;

public:
    /// Constructor
    IsotropicDamageMaterial(int n, Domain *d);
    /// Destructor
    virtual ~IsotropicDamageMaterial();

    bool hasMaterialModeCapability(MaterialMode mode) const override;
    const char *giveClassName() const override { return "IsotropicDamageMaterial"; }

    /// Returns reference to undamaged (bulk) material
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep) override;

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<IsotropicDamageMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<4> giveRealStressVector_PlaneStrain( const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<IsotropicDamageMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArray giveRealStressVector_StressControl(const FloatArray &strain, const IntArray &strainControl, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<IsotropicDamageMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<IsotropicDamageMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<IsotropicDamageMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    FloatArrayF<6> giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override;
    virtual double evaluatePermanentStrain(double kappa, double omega) const { return 0.; }

    /**
     * Returns the value of material property 'aProperty'. Property must be identified
     * by unique int id. Integration point also passed to allow for materials with spatially
     * varying properties
     * @param aProperty ID of property requested.
     * @param gp Integration point,
     * @return Property value.
     */
    double give(int aProperty, GaussPoint *gp) const override;
    /**
     * Computes the equivalent strain measure from given strain vector (full form).
     * @param[out] kappa Return parameter, containing the corresponding equivalent strain.
     * @param strain Total strain vector in full form.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual double computeEquivalentStrain(const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) const = 0;
    /**Computes derivative of the equivalent strain with regards to strain
     * @param[out] answer Contains the resulting derivative.
     * @param strain Strain vector.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void computeEta(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) const { OOFEM_ERROR("not implemented"); }
    /**
     * Computes the value of damage parameter omega, based on given value of equivalent strain.
     * @param[out] omega Contains result.
     * @param kappa Equivalent strain measure.
     * @param strain Total strain in full form.
     * @param gp Integration point.
     */
    virtual double computeDamageParam(double kappa, const FloatArray &strain, GaussPoint *gp) const = 0;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new IsotropicDamageMaterialStatus(gp); }

    FloatMatrixF<1,1> give1dStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp,
                                            TimeStep *tStep) const override;
  
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
    
protected:
    /**
     * Abstract service allowing to perform some initialization, when damage first appear.
     * @param kappa Scalar measure of strain level.
     * @param totalStrainVector Current total strain vector.
     * @param gp Integration point.
     */
    virtual void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp) const { }

    /**
     * Returns the value of derivative of damage function
     * wrt damage-driving variable kappa corresponding
     * to a given value of the  kappa, depending on
     * the type of selected damage law.
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    virtual double damageFunctionPrime(double kappa, GaussPoint *gp) const {
        OOFEM_ERROR("not implemented");
        return 0;
    }

    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp,TimeStep *tStep) const override;
    FloatMatrixF<4,4> givePlaneStrainStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;

};
} // end namespace oofem
#endif // isodamagemodel_h
