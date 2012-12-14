
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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef latticedamage2d_h

#include "material.h"
#include "linearelasticmaterial.h"
#include "latticematstatus.h"
#include "cltypes.h"
#include "structuralmaterial.h"
#include "randommaterialext.h"
#include "strainvector.h"
#include "stressvector.h"

namespace oofem {
/**
 * This class implements associated Material Status to LatticeDamage2d.
 */
class LatticeDamage2dStatus : public LatticeMaterialStatus, public RandomMaterialStatusExtensionInterface
{
protected:
    /// scalar measure of the largest strain level ever reached in material
    double kappa;
    /// non-equilibrated scalar measure of the largest strain level
    double tempKappa;

    /// scalar measure of the strain
    double equivStrain;
    /// non-equilibrated scalar measure of the strain
    double tempEquivStrain;

    /// damage level of material
    double damage;
    /// non-equilibrated damage level of material
    double tempDamage;

    /// dissipation
    double dissipation;
    /// non-equilibrated dissipation
    double tempDissipation;

    /// increment of dissipation
    double deltaDissipation;
    /// non-equilibrated increment of dissipation
    double tempDeltaDissipation;

    /// characteristic length
    double le;

    /// random material parameter stored in status, since each gp has a differnet value.
    double e0;

    /** the crack_flag indicates if the gp is cracked:
     *  crack_flag = 0 gp is uncracked
     * crack_flag = 1 gp is cracked and damage grows
     * crack_flag = 2 gp is cracked and damage does not grow
     */
    int crack_flag;

    /// non-equilibrated temp flag
    int temp_crack_flag;

    /// crack width
    double crackWidth;

    /// non-equilibrated crack width
    double tempCrackWidth;

    /// reduced strain
    FloatArray reducedStrain;

    /// non-equilibrated reduced strain
    FloatArray tempReducedStrain;

public:

    /// Constructor
    LatticeDamage2dStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~LatticeDamage2dStatus() {}


    /// Returns the last equilibrated scalar measure of the largest strain level
    double giveKappa() { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level
    double giveTempKappa() { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value
    void setTempKappa(double newKappa) { tempKappa = newKappa; }

    /// Returns the last equilibrated scalar measure of the largest strain level
    double giveEquivalentStrain() { return equivStrain; }
    /// Returns the temp. scalar measure of the largest strain level
    double giveTempEquivalentStrain() { return tempEquivStrain; }
    /// Sets the temp scalar measure of the largest strain level to given value
    void setTempEquivalentStrain(double newEquivStrain) { tempEquivStrain = newEquivStrain; }

    /// Returns the last equilibrated dissipation
    double giveDissipation() { return dissipation; }
    /// Returns the temp. dissipation
    double giveTempDissipation() { return tempDissipation; }
    /// Sets the temp dissipation
    void setTempDissipation(double newDiss) { tempDissipation = newDiss; }

    /// Returns the last equilibrated increment of dissipation
    double giveDeltaDissipation() { return deltaDissipation; }
    /// Returns the temp. increment dissipation
    double giveTempDeltaDissipation() { return tempDeltaDissipation; }
    /// Sets the temp. increment dissipation
    void setTempDeltaDissipation(double newDiss) { tempDeltaDissipation = newDiss; }


    /// Returns the last equilibrated damage level
    double giveDamage() { return damage; }
    /// Returns the temp. damage level
    double giveTempDamage() { return tempDamage; }
    /// Sets the temp damage level to given value
    void setTempDamage(double newDamage) { tempDamage = newDamage; }

    /*
     * Assign the temp value of plastic strain.
     * @v new temp value of plastic strain
     */
    void giveTempReducedStrain(FloatArray &answer) const
    { answer = tempReducedStrain; }

    void giveReducedStrain(FloatArray &answer) const
    { answer = reducedStrain; }

    /*
     * Assign the temp value of plastic strain.
     * @v new temp value of plastic strain
     */
    void letTempReducedStrainBe(const FloatArray &v)
    { tempReducedStrain = v; }

    /// Prints the receiver state to given stream
    void printOutputAt(FILE *file, TimeStep *tStep);

    /// Returns characteristic length stored in receiver
    double giveLe() { return le; }

    /// Sets characteristic length to given value
    void setLe(double ls) { le = ls; }

    ///Sets the temp_crack_flag
    void setTempCrackFlag(int val);

    /// Returns the crack_flag
    int giveCrackFlag();

    /// Set random e0
    void setE0(double val) { e0 = val; }

    ///Sets the temp crack width
    void setTempCrackWidth(double val);
    //Gives the last equilibrated crack width
    double giveCrackWidth() { return this->crackWidth; }

    // definition
    virtual const char *giveClassName() const { return "LatticeDamage2dStatus"; }
    virtual classType giveClassID() const { return LatticeDamage2dStatusClass; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached

    virtual Interface *giveInterface(InterfaceType);

    virtual void setVariableInStatus(double variable);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};


/**
 * This class implements a local random isotropic damage model for concrete in tension for 2D lattice elements.
 */
class LatticeDamage2d : public StructuralMaterial, public RandomMaterialExtensionInterface
    //
{
protected:

    /// Normal modulus
    double eNormal;
    /// Shear modulus
    double eShear;
    /// Torsion modulus
    double eTorsion;
    /// ratio of shear and normal modulus
    double alphaOne;
    /// ratio of torsion and normal modulus
    double alphaTwo;

    /// mean effective strain at peak
    double e0Mean;

    /// mean effective strain at sigma1
    double e0OneMean;

    /**parameter which determines the typ of the softeningFunction
     * 1 = linear softening
     * 2 = bilinear softening
     * 3 = exponential softening
     **/
    int softeningType;

    /// determines the softening -> corresponds to threshold of crack opening (not strain)
    double wf, wfOne;

    /// parameter for the elliptic equivalent strain function
    double ec;

    /// parameter setting ratio of shear and tensile strength
    double coh;

    /// coefficient variation of the Gaussian distribution
    double coefficientOfVariation;

    /// flag which chooses between no distribution (0) and Gaussian distribution (1)
    double localRandomType;

public:

    /// Constructor
    LatticeDamage2d(int n, Domain *d);
    /// Destructor
    virtual ~LatticeDamage2d();

    virtual const char *giveClassName() const { return "LatticeDamage2d"; }

    virtual classType giveClassID() const { return LatticeDamage2dClass; }

    virtual int giveStressStrainComponentIndOf(MatResponseForm, MaterialMode mmode, int);

    virtual void giveStressStrainMask(IntArray & answer, MatResponseForm, MaterialMode mmode) const;

    virtual int giveSizeOfReducedStressStrainVector(MaterialMode);

    virtual void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *,
                                                 const FloatArray &charVector3d);

    virtual void giveFullCharacteristicVector(FloatArray &answer,  GaussPoint *,
                                              const FloatArray &);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    /**
     * Computes the tangent stiffness.
     * @param answer return param, containing the stiffness
     * @param gp integration point
     * @param atTime timeStep
     */

    void giveTangentStiffnessMatrix(FloatMatrix &answer,
                                    GaussPoint *gp,
                                    TimeStep *atTime);

    /**
     * Computes the secant stiffness.
     * @param answer return param, containing the stiffness
     * @param gp integration point
     * @param atTime timeStep
     */

    void giveSecantStiffnessMatrix(FloatMatrix &answer,
                                   GaussPoint *gp,
                                   TimeStep *atTime);


    /**
     * Computes the elastic stiffness.
     * @param answer return param, containing the stiffness
     * @param gp integration point
     * @param atTime timeStep
     */

    void giveElasticStiffnessMatrix(FloatMatrix &answer,
                                    GaussPoint *gp,
                                    TimeStep *atTime);

    virtual int hasMaterialModeCapability(MaterialMode mode);

    /**
     * Computes the equivalent strain measure from given strain vector (full form).
     * @param kappa return param, comtaining the corresponding equivalent strain
     * @param strain total strain vector in full form
     * @param gp integration point
     * @param atTime timeStep
     */
    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime);

    /**
     * computes the value of damage parameter omega, based on given value of equivalent strain
     * @param omega contains result
     * @param kappa equivalent strain measure
     */
    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp);
    /** Interface requesting service */
    virtual Interface *giveInterface(InterfaceType);


    virtual void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                                      const FloatArray &, TimeStep *);

    /** Reimplemented from RandomMaterialInterface */
    virtual void giveRandomParameters(FloatArray &param);

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp,  TimeStep *tStep);

    /// Creates corresponding status
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;


    virtual MaterialStatus *giveStatus(GaussPoint *gp) const;

    virtual double give(int aProperty, GaussPoint *gp);


protected:


    /**
     *  Performs initialization, when damage first appear. The Le characteristic length is
     *  computed from the direction of largest positive principal strain and stored
     *  in corresponding status.
     *  @param kappa scalar measure of strain level
     *  @param totalStrainVector current total strain vector
     *  @param gp integration point
     */
    void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp);


    virtual int giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *atTime);

    virtual int giveIPValueSize(InternalStateType type,
                                GaussPoint *gp);

    virtual int giveIntVarCompFullIndx(IntArray &answer,
                                       InternalStateType type,
                                       MaterialMode mmode);

    virtual InternalStateValueType giveIPValueType(InternalStateType type);
};
} // end namespace oofem
#define latticedamage2d_h
#endif
