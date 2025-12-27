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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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
    double kappa = 0.;
    /// Non-equilibrated scalar measure of the largest equivalent displacement.
    double tempKappa = 0.;
    /// Damage level of material.
    double damage = 0.;
    /// Non-equilibrated damage level of material.
    double tempDamage = 0.;
public:
    /// Constructor
    IntMatIsoDamageStatus(GaussPoint *g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    /// Returns the last equilibrated scalar measure of the largest jump level.
    double giveKappa() const { return kappa; }
    /// Returns the temp. scalar measure of the largest jump level.
    double giveTempKappa() const { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value.
    void setTempKappa(double newKappa) { tempKappa = newKappa; }
    /// Returns the last equilibrated damage level.
    double giveDamage() const override { return damage; }
    /// Returns the temp. damage level.
    double giveTempDamage() const override { return tempDamage; }
    /// Sets the temp damage level to given value.
    void setTempDamage(double newDamage) { tempDamage = newDamage; }

    const char *giveClassName() const override { return "IntMatIsoDamageStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
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
    double kn = 0.;
    /// Shear moduli.
    double ks = 0.;
    /// Tension strength.
    double ft = 0.;
    /// Fracture energy.
    double gf = 0.;
    /// Limit elastic deformation.
    double e0 = 0.;
    /// Maximum limit on omega. The purpose is elimination of a too compliant material which may cause convergency problems. Set to something like 0.99 if needed.
    double maxOmega = 0.999999;

    bool semiExplicit = false; // If semi-explicit time integration should be used

public:
    IntMatIsoDamage(int n, Domain *d);

    const char *giveInputRecordName() const override { return _IFT_IntMatIsoDamage_Name; }
    const char *giveClassName() const override { return "IntMatIsoDamage"; }

    bool hasAnalyticalTangentStiffness() const override { return true; } 

    FloatArrayF<3> giveEngTraction_3d(const FloatArrayF<3> &jump, GaussPoint *gp,TimeStep *tStep) const override;

    FloatArrayF<3> giveFirstPKTraction_3d(const FloatArrayF<3> &jump, const FloatMatrixF<3,3> &F, GaussPoint *gp, TimeStep *tStep) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    /**
     * Computes the equivalent jump measure from given jump vector (full form).
     * @param jump Jump vector in full form.
     * @param gp Integration point.
     * @param tStep Time step.
     * @return Return parameter containing the corresponding equivalent jump (kappa)
     */
    virtual double computeEquivalentJump(const FloatArray &jump) const;

    /**
     * computes the value of damage parameter omega, based on given value of equivalent strain.
     * @param kappa Equivalent strain measure.
     * @param strain Total strain vector in full form. (unnecessary?)
     * @param gp Integration point.
     * @return omega.
     */
    virtual double computeDamageParam(double kappa) const;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<IntMatIsoDamageStatus>(gp); }

    FloatMatrixF<2,2> give2dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> give3dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;
};
} // end namespace oofem
#endif // isointerfacedamage01_h
