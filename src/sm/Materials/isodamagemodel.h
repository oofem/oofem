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
#include "Materials/linearelasticmaterial.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"

///@name Input fields for IsotropicDamageMaterial
//@{
#define _IFT_IsotropicDamageMaterial_talpha "talpha"
#define _IFT_IsotropicDamageMaterial_maxOmega "maxomega"
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
    double kappa;
    /// Non-equilibrated scalar measure of the largest strain level.
    double tempKappa;
    /// Damage level of material.
    double damage;
    /// Non-equilibrated damage level of material.
    double tempDamage;
    /**
     * Characteristic element length,
     * computed when damage initialized from direction of
     * maximum positive principal strain. Fixed during further loading.
     */
    double le;
    /// Angle characterizing the crack direction.
    double crack_angle;
    /// Crack orientation normalized to damage magnitude. This is useful for plotting cracks as a vector field (paraview etc.).
    FloatArray crackVector;

#ifdef keep_track_of_dissipated_energy
    /// Density of total work done by stresses on strain increments.
    double stressWork;
    /// Non-equilibrated density of total work done by stresses on strain increments.
    double tempStressWork;
    /// Density of dissipated work.
    double dissWork;
    /// Non-equilibrated density of dissipated work.
    double tempDissWork;
#endif

public:
    /// Constructor
    IsotropicDamageMaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~IsotropicDamageMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    /// Returns the last equilibrated scalar measure of the largest strain level.
    double giveKappa() { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level.
    double giveTempKappa() { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value.
    void setTempKappa(double newKappa) { tempKappa = newKappa; }
    /// Returns the last equilibrated damage level.
    double giveDamage() { return damage; }
    /// Returns the temp. damage level.
    double giveTempDamage() { return tempDamage; }
    /// Sets the temp damage level to given value.
    void setTempDamage(double newDamage) { tempDamage = newDamage; }

    /// Returns characteristic length stored in receiver.
    double giveLe() { return le; }
    /// Sets characteristic length to given value.
    void setLe(double ls) { le = ls; }
    /// Returns crack angle stored in receiver.
    double giveCrackAngle() { return crack_angle; }
    /// Sets crack angle to given value.
    void setCrackAngle(double ca) { crack_angle = ca; }
    /// Returns crack vector stored in receiver. This is useful for plotting cracks as a vector field (paraview etc.).
    void giveCrackVector(FloatArray &answer);
    /// Sets crack vector to given value. This is useful for plotting cracks as a vector field (paraview etc.).
    void setCrackVector(FloatArray cv) { crackVector = cv; }

#ifdef keep_track_of_dissipated_energy
    /// Returns the density of total work of stress on strain increments.
    double giveStressWork() { return stressWork; }
    /// Returns the temp density of total work of stress on strain increments.
    double giveTempStressWork() { return tempStressWork; }
    /// Sets the density of total work of stress on strain increments to given value.
    void setTempStressWork(double w) { tempStressWork = w; }
    /// Returns the density of dissipated work.
    double giveDissWork() { return dissWork; }
    /// Returns the density of temp dissipated work.
    double giveTempDissWork() { return tempDissWork; }
    /// Sets the density of dissipated work to given value.
    void setTempDissWork(double w) { tempDissWork = w; }
    /// Computes the increment of total stress work and of dissipated work.
    void computeWork(GaussPoint *gp);
#endif

    // definition
    virtual const char *giveClassName() const { return "IsotropicDamageMaterialModelStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
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
    double tempDillatCoeff;

    /// Maximum limit on omega. The purpose is elimination of a too compliant material which may cause convergence problems. Set to something like 0.99 if needed.
    double maxOmega;

    /// Reference to bulk (undamaged) material
    LinearElasticMaterial *linearElasticMaterial;
    /**
     * Variable controlling type of loading/unloading law, default set to idm_strainLevel
     * defines the two two possibilities:
     * - idm_strainLevelCR the unloaing takes place, when strain level is smaller than the largest level ever reached;
     * - idm_damageLevelCR the unloaing takes place, when damage level is smaller than the largest damage ever  reached;
     */
    enum loaUnloCriterium { idm_strainLevelCR, idm_damageLevelCR } llcriteria;

public:
    /// Constructor
    IsotropicDamageMaterial(int n, Domain * d);
    /// Destructor
    virtual ~IsotropicDamageMaterial();

    virtual int hasNonLinearBehaviour() { return 1; }

    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual const char *giveClassName() const { return "IsotropicDamageMaterial"; }

    /// Returns reference to undamaged (bulk) material
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer,  GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_StressControl(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, const IntArray &strainControl, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *, TimeStep *);
    /**
     * Returns the value of material property 'aProperty'. Property must be identified
     * by unique int id. Integration point also passed to allow for materials with spatially
     * varying properties
     * @param aProperty ID of property requested.
     * @param gp Integration point,
     * @return Property value.
     */
    virtual double give(int aProperty, GaussPoint *gp);
    /**
     * Computes the equivalent strain measure from given strain vector (full form).
     * @param[out] kappa Return parameter, containing the corresponding equivalent strain.
     * @param strain Total strain vector in full form.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) = 0;
    /**Computes derivative of the equivalent strain with regards to strain
     * @param[out] answer Contains the resulting derivative.
     * @param strain Strain vector.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void computeEta(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) { OOFEM_ERROR("not implemented"); }
    /**
     * Computes the value of damage parameter omega, based on given value of equivalent strain.
     * @param[out] omega Contains result.
     * @param kappa Equivalent strain measure.
     * @param strain Total strain in full form.
     * @param gp Integration point.
     */
    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp) = 0;

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new IsotropicDamageMaterialStatus(1, domain, gp); }

protected:
    /**
     * Abstract service allowing to perform some initialization, when damage first appear.
     * @param kappa Scalar measure of strain level.
     * @param totalStrainVector Current total strain vector.
     * @param gp Integration point.
     */
    virtual void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp) { }

    /**
     * Returns the value of derivative of damage function
     * wrt damage-driving variable kappa corresponding
     * to a given value of the  kappa, depending on
     * the type of selected damage law.
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    virtual double damageFunctionPrime(double kappa, GaussPoint *gp) {
        OOFEM_ERROR("not implemented");
        return 0;
    }

    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mmode,
                                       GaussPoint *gp,
                                       TimeStep *tStep);
};
} // end namespace oofem
#endif // isodamagemodel_h
