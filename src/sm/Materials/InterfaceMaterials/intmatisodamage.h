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

#ifndef intmatisodamage_h
#define intmatisodamage_h

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"


///@name Input fields for IntMatIsoDamage
//@{
#define _IFT_IntMatIsoDamage_Name "intmatisodamage"
#define _IFT_IntMatIsoDamage_kn "kn" /// Normal stiffness
#define _IFT_IntMatIsoDamage_ks "ks" /// Tangent stiffness
#define _IFT_IntMatIsoDamage_ft "ft"
#define _IFT_IntMatIsoDamage_gf "gf"
#define _IFT_IntMatIsoDamage_maxOmega "maxomega"

//@}

namespace oofem {
class GaussPoint;

/**
 * This class implements the InterfaceMaterialStatus associated with IntMatIsoDamage.
 */
class IntMatIsoDamageStatus : public StructuralInterfaceMaterialStatus
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
    IntMatIsoDamageStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~IntMatIsoDamageStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    /// Returns the last equilibrated scalar measure of the largest jump level.
    double giveKappa() { return kappa; }
    /// Returns the temp. scalar measure of the largest jump level.
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
    virtual const char *giveClassName() const { return "IntMatIsoDamageStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};


/**
 * Simple isotropic damage based model for interface elements.
 * In 2d, the interface elements are used to model contact layer between
 * element edges. The jump vector contains the relative displacements
 * (in normal and shear direction). The traction vector contains the corresponding
 * tractions in normal and tangent direction.
 *
 * The behaviour of the model is elastic, described by normal and shear stiffness components.
 * Isotropic damage is initiated when the stress reaches the tensile strength. The damage evolution
 * is governed by the normal jump
 */
class IntMatIsoDamage : public StructuralInterfaceMaterial
{
protected:
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
    /// Maximum limit on omega. The purpose is elimination of a too compliant material which may cause convergency problems. Set to something like 0.99 if needed.
    double maxOmega;

    bool semiExplicit; // If semi-explicit time integration should be used

public:
    /// Constructor
    IntMatIsoDamage(int n, Domain *d);
    /// Destructor
    virtual ~IntMatIsoDamage();

    virtual const char *giveInputRecordName() const { return _IFT_IntMatIsoDamage_Name; }
    virtual bool hasAnalyticalTangentStiffness() const { return true; } 

    virtual void giveEngTraction_3d(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &jump, TimeStep *tStep);

    virtual void giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                        const FloatMatrix &F, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    

    /**
     * Computes the equivalent jump measure from given jump vector (full form).
     * @param[out] kappa Return parameter containing the corresponding equivalent jump.
     * @param jump Jump vector in full form.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void computeEquivalentJump(double &kappa, const FloatArray &jump);

    /**
     * computes the value of damage parameter omega, based on given value of equivalent strain.
     * @param[out] omega Contains result.
     * @param kappa Equivalent strain measure.
     * @param strain Total strain vector in full form. (unnecessary?)
     * @param gp Integration point.
     */
    virtual void computeDamageParam(double &omega, double kappa);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new IntMatIsoDamageStatus(1, domain, gp); }

protected:
    virtual void give2dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode,
                                                GaussPoint *gp, TimeStep *tStep);
    virtual void give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode,
                                                GaussPoint *gp, TimeStep *tStep);
};
} // end namespace oofem
#endif // isointerfacedamage01_h
