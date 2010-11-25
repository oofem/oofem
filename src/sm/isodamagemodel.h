/* $Header: /home/cvs/bp/oofem/sm/src/isodamagemodel.h,v 1.5 2003/04/06 14:08:30 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   ******************************************************
//   *** CLASS SIMPLE ISOTROPIC DAMAGE MODEL   ************
//   ******************************************************

#ifndef isodamagemodel_h
#define isodamagemodel_h

// this turns on or off a bunch of internal variables
// that allow tracing the distribution of dissipated energy
// (can be turned off if such information is not needed)
#define keep_track_of_dissipated_energy

#include "material.h"
#include "linearelasticmaterial.h"
#include "structuralmaterial.h"
#include "structuralms.h"

namespace oofem {
// material contant's keys for give()
class GaussPoint;

/**
 * This class implements associated Material Status to IsotropicDamageMaterial.
 * It is atribute of matStatusDictionary at every GaussPoint, for which this material
 * is active.
 */
class IsotropicDamageMaterialStatus : public StructuralMaterialStatus
{
    /*
     * This class implements associated Material Status to IsotropicDamageMaterial.
     * It is atribute of matStatusDictionary at every GaussPoint, for which this material
     * is active.
     * DESCRIPTION:
     * Idea used there is that we have variables
     * describing:
     *   1) state at previous equilibrium state (variables without temp)
     * 2) state during searching new equilibrium (variables with temp)
     * when we start search new state from previous equilibrium one we copy
     * non-tem variables into temp ones. And after we reach new equilibrium
     * (now decribed by temp variables) we copy tem-var into non-tepm ones
     * (see function updateYourself).
     *
     * variables description:
     *
     * kappa      - scalar measure of the largest strain level ever reached in material
     * tempKappa
     *
     */

protected:
    /// scalar measure of the largest strain level ever reached in material
    double kappa;
    /// non-equilibrated scalar measure of the largest strain level
    double tempKappa;
    /// damage level of material
    double damage;
    /// non-equilibrated damage level of material
    double tempDamage;

#ifdef keep_track_of_dissipated_energy
    /// density of total work done by stresses on strain increments
    double stressWork;
    /// non-equilibrated density of total work done by stresses on strain increments
    double tempStressWork;
    /// density of dissipated work
    double dissWork;
    /// non-equilibrated density of dissipated work
    double tempDissWork;
#endif

public:
    /// Constructor
    IsotropicDamageMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    ~IsotropicDamageMaterialStatus();

    /// Prints the receiver state to stream
    void   printOutputAt(FILE *file, TimeStep *tStep);

    /// Returns the last equilibrated scalar measure of the largest strain level
    double giveKappa() { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level
    double giveTempKappa() { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value
    void   setTempKappa(double newKappa) { tempKappa = newKappa; }
    /// Returns the last equilibrated damage level
    double giveDamage() { return damage; }
    /// Returns the temp. damage level
    double giveTempDamage() { return tempDamage; }
    /// Sets the temp damage level to given value
    void   setTempDamage(double newDamage) { tempDamage = newDamage; }

#ifdef keep_track_of_dissipated_energy
    /// Returns the density of total work of stress on strain increments
    double giveStressWork() { return stressWork; }
    /// Returns the temp density of total work of stress on strain increments
    double giveTempStressWork() { return tempStressWork; }
    /// Sets the density of total work of stress on strain increments to given value
    void setTempStressWork(double w) { tempStressWork = w; }
    /// Returns the density of dissipated work
    double giveDissWork() { return dissWork; }
    /// Returns the density of temp dissipated work
    double giveTempDissWork() { return tempDissWork; }
    /// Sets the density of dissipated work to given value
    void setTempDissWork(double w) { tempDissWork = w; }
    /// computes the increment of total stress work and of dissipated work
    void computeWork(GaussPoint *);
#endif

    // definition
    const char *giveClassName() const { return "IsotropicDamageMaterialModelStatus"; }
    classType             giveClassID() const { return IsotropicDamageMaterialStatusClass; }

    /**
     * Initializes the temporary internal variables, describing the current state according to
     * previously reached equilibrium internal variables.
     */
    virtual void initTempStatus();
    /**
     * Update equilibrium history variables according to temp-variables.
     * Invoked, after new equilibrium state has been reached.
     */
    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached

    // saves current context(state) into stream
    /**
     * Stores context of receiver into given stream.
     * Only non-temp internal history variables are stored.
     * @param stream stream where to write data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores context of receiver from given stream.
     * @param stream stream where to read data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};


/**
 * Base class representing general isotropic damage model.
 * It is based on isotropic damage concept, assuming that damage evolution law
 * is postulated in explicit form, relatin damage parameter (omega) to scalar measure
 * of the largest strain level ever reached in material (kappa).
 */
class IsotropicDamageMaterial : public StructuralMaterial
{
    /*
     *
     * DESCRIPTION
     * This class implements a general isotropic damage model
     *
     * A model is based on isotropic damage concept, assuming that damage evolution law
     * is postulated in explicit form, relatin damage parameter (omega) to scalar measure
     * of the largest strain level ever reached in material (kappa).
     *
     * TASK
     * - Returning standard material stiffness and flexibility marices for 3d-case.
     * according to current state determined by using data stored
     * in Gausspoint.
     * - Returning a material property (method 'give'). Only for non-standard elements.
     * - Returning real stress state vector(tensor) at gauss point for 3d - case.
     * - Storing & restoring Material Status sored in gp matStatusDictionary.
     */

protected:
    /// Coefficient of thermal dilatation
    double tempDillatCoeff;

    /// Maximum limit on omega. The purpose is elimination of a too compliant material which may cause convergency problems. Set to something like 0.99 if needed
    double maxOmega;

    /// Reference to bulk (undamaged) material
    LinearElasticMaterial *linearElasticMaterial;
    /**
     * variable controling type of loading/unloading law, default set to idm_strainLevel
     * defines the two two possibilities:
     * - idm_strainLevelCR the unloaing takes place, when strain level is smaller than the largest level ever reached;
     * - idm_damageLevelCR the unloaing takes place, when damage level is smaller than the largest damage ever  reached;
     */
    enum loaUnloCriterium { idm_strainLevelCR, idm_damageLevelCR } llcriteria;

public:
    /// Constructor
    IsotropicDamageMaterial(int n, Domain *d);
    /// Destructor
    ~IsotropicDamageMaterial();

    /// Returns nonzero indicating that receiver is nonlinear
    int hasNonLinearBehaviour()   { return 1; }
    /**
     * Tests, if material supports material mode.
     * @param mode required material mode
     * @return nonzero if supported, zero otherwise
     */
    int hasMaterialModeCapability(MaterialMode mode);
    const char *giveClassName() const { return "IsotropicDamageMaterial"; }
    classType giveClassID()         const { return IsotropicDamageMaterialClass; }

    /// Returns reference to undamaged (bulk) material
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }
    /**
     * Computes full 3d material stiffness matrix at given integration point, time, respecting load history
     * in integration point.
     * @param answer computed results
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                               MatResponseForm, MatResponseMode,
                                               GaussPoint * gp,
                                               TimeStep * atTime);

    /**
     * Computes the real stress vector for given strain increment and integration point.
     * Temp vars are updated accordingly
     * @param answer contains result
     * @param form material response form
     * @param gp integration point
     * @param reducedStrain strain  vector in reduced form
     * @param tStep current time step (most models are able to respond only when atTime is current time step)
     */
    void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);


#ifdef __OOFEG
#endif
    /**
     * Returns the integration point corresponding value in Reduced form.
     * @param answer contain corresponding ip value, zero sized if not available
     * @param aGaussPoint integration point
     * @param type determines the type of internal variable
     * @param type determines the type of internal variable
     * @returns nonzero if ok, zero if var not supported
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
     * @param type determines the internal variable requested (physical meaning)
     * @returns nonzero if ok or error is generated for unknown mat mode.
     */
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);


    /**
     * Returns the type of internal variable (scalar, vector, tensor,...).
     * @param type determines the type of internal variable
     * @returns type of internal variable
     */
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    /**
     * Returns the corresponding integration point  value size in Reduced form.
     * @param type determines the type of internal variable
     * @returns var size, zero if var not supported
     */
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
    /**
     * Returns a vector of coefficients of thermal dilatation in direction
     * of each material principal (local) axis.
     * @param answer vector of thermal dilatation coefficients
     * @param gp integration point
     * @param tStep time step (most models are able to respond only when atTime is current time step)
     */
    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *, TimeStep *);



    /**
     * Computes the equivalent strain measure from given strain vector (full form).
     * @param kappa return param, comtaining the corresponding equivalent strain
     * @param strain total strain vector in full form
     * @param gp integration point
     * @param atTime timestep
     */
    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime) = 0;
    /**
     * computes the value of damage parameter omega, based on given value of equivalent strain
     * @param omega contains result
     * @param kappa equivalent strain measure
     */
    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp) = 0;
    /**
     * Instanciates the receiver from input record
     */
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    /// Creates corresponding material status
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new IsotropicDamageMaterialStatus(1, domain, gp); }

protected:

    /**
     * Abstract service allowing to perfom some initialization, when damage first appear
     * @param kappa scalar measure of strain level
     * @param totalStrainVector current total strain vector
     * @param gp integration point
     */
    virtual void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp) { }


    // Overloaded to use specialized versions of these services possibly implemented by linearElastic member
    /**
     * Method for computing plane stress stifness matrix of receiver.
     * Default implementation overloaded to use direct implementation of
     * corresponding service at bulk material model level.
     * @param answer stifness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    void givePlaneStressStiffMtrx(FloatMatrix & answer, MatResponseForm form, MatResponseMode,
                                  GaussPoint * gp,
                                  TimeStep * atTime);
    /**
     * Method for computing plane strain stifness matrix of receiver.
     * Default implementation overloaded to use direct implementation of
     * corresponding service at bulk material model level.
     * @param answer stifness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    void givePlaneStrainStiffMtrx(FloatMatrix & answer, MatResponseForm form, MatResponseMode,
                                  GaussPoint * gp,
                                  TimeStep * atTime);
    /**
     * Method for computing 1d stifness matrix of receiver.
     * Default implementation overloaded to use direct implementation of
     * corresponding service at bulk material model level.
     * @param answer stifness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    void give1dStressStiffMtrx(FloatMatrix & answer, MatResponseForm form, MatResponseMode,
                               GaussPoint * gp,
                               TimeStep * atTime);
};
} // end namespace oofem
#endif // isodamagemodel_h
