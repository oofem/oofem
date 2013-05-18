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

#ifndef slipmodelcebfip90_h
#define slipmodelcebfip90_h

#include "material.h"
#include "linearelasticmaterial.h"
#include "structuralmaterial.h"
#include "structuralms.h"

namespace oofem {

/**
 * This class implements associated Material Status to IsoInterfaceDamageMaterial.
 */
class IsoInterfaceDamageMaterialStatus : public StructuralMaterialStatus
{
protected:
    /// Scalar measure of the largest equivalent displacement ever reached in material.
    double kappa;
    /// Non-equilibrated scalar measure of the largest equivalent displacement.
    double tempKappa;
    /// Damage level of material.
    double damage;
    /// Non-equilibrated damage level of material.
    double tempDamage;

public:
    /// Constructor
    IsoInterfaceDamageMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~IsoInterfaceDamageMaterialStatus();

    void printOutputAt(FILE *file, TimeStep *tStep);

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

    // definition
    virtual const char *giveClassName() const { return "IsoInterfaceDamageMaterialStatus"; }
    virtual classType giveClassID() const { return MaterialStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    // saves current context(state) into stream
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};


/**
 * Base class representing general isotropic damage model for interfaces.
 * It is based on isotropic damage concept, assuming that damage evolution law
 * is postulated in explicit form, relation damage parameter (omega) to scalar measure
 * of the largest strain level ever reached in material (kappa).
 */
class IsoInterfaceDamageMaterial : public StructuralMaterial
{
protected:
    /// Coefficient of thermal dilatation.
    double tempDillatCoeff;
    /// Elastic properties (normal moduli).
    double kn;
    /// Shear moduli.
    double ks;
    /// Tension strength.
    double ft;
    /// Fracture energy.
    double gf;
    /// Limit elastic deformation.
    double e0;

public:
    /// Constructor
    IsoInterfaceDamageMaterial(int n, Domain *d);
    /// Destructor
    virtual ~IsoInterfaceDamageMaterial();

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual const char *giveClassName() const { return "IsoInterfaceDamageMaterial"; }
    virtual classType giveClassID() const { return StructuralMaterialClass; }

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer,  MatResponseForm form, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep);

    virtual int giveStressStrainComponentIndOf(MatResponseForm form, MaterialMode mmode, int ind);
    virtual void giveStressStrainMask(IntArray &answer, MatResponseForm form, MaterialMode mmode) const;
    virtual int giveSizeOfReducedStressStrainVector(MaterialMode mmode);
    void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                         const FloatArray &charVector3d);
    void giveFullCharacteristicVector(FloatArray &answer,  GaussPoint *gp,
                                      const FloatArray &strainVector);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);

    virtual InternalStateValueType giveIPValueType(InternalStateType type);

    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new IsoInterfaceDamageMaterialStatus(1, domain, gp); }

protected:
    // Overloaded to use specialized versions of these services possibly implemented by linearElastic member
    void give2dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode,
                                                GaussPoint *gp, TimeStep *tStep);
};
} // end namespace oofem
#endif // slipmodelcebfip90_h
