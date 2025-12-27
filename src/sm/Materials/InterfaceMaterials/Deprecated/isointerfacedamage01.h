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

#ifndef isointerfacedamage01_h
#define isointerfacedamage01_h

#include "sm/Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"

///@name Input fields for IsoInterfaceDamageMaterial
//@{
#define _IFT_IsoInterfaceDamageMaterial_Name "isointrfdm01"
#define _IFT_IsoInterfaceDamageMaterial_kn "kn"
#define _IFT_IsoInterfaceDamageMaterial_ks "ks"
#define _IFT_IsoInterfaceDamageMaterial_ft "ft"
#define _IFT_IsoInterfaceDamageMaterial_gf "gf"
#define _IFT_IsoInterfaceDamageMaterial_maxOmega "maxomega"
#define _IFT_IsoInterfaceDamageMaterial_beta "beta"
//@}

namespace oofem {

/**
 * This class implements associated Material Status to IsoInterfaceDamageMaterial.
 */
class IsoInterfaceDamageMaterialStatus : public StructuralInterfaceMaterialStatus
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
    IsoInterfaceDamageMaterialStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    /// Returns the last equilibrated scalar measure of the largest strain level.
    double giveKappa() const { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level.
    double giveTempKappa() const { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value.
    void setTempKappa(double newKappa) { tempKappa = newKappa; }
    /// Returns the last equilibrated damage level.
    double giveDamage() const override { return damage; }
    /// Returns the temp. damage level.
    double giveTempDamage() const override { return tempDamage; }
    /// Sets the temp damage level to given value.
    void setTempDamage(double newDamage) { tempDamage = newDamage; }

    const char *giveClassName() const override { return "IsoInterfaceDamageMaterialStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * Simple isotropic damage based model for 2d interface elements.
 * In 2d, the interface elements are used to model contact layer between
 * element edges. The generalized strain vector contains two relative displacements
 * (in normal and shear direction). The generalized stress vector contains corresponding
 * tractions in normal and tangent direction.
 *
 * The behaviour of the model is elastic, described by normal and shear stiffness components.
 * Isotropic damage is initiated  when the stress reaches the tensile strength. Damage evolution
 * is governed by normal component of generalized strain vector (normal relative displacement)
 * by an exponential softening law.
 */
class IsoInterfaceDamageMaterial : public StructuralInterfaceMaterial
{
protected:
    /// Coefficient of thermal dilatation.
    double tempDillatCoeff = 0.;
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
    /// Weight factor for the influence of shear component of displacement jump on equivalent strain.
    double beta = 0.;

public:
    /// Constructor
    IsoInterfaceDamageMaterial(int n, Domain * d);

    bool hasAnalyticalTangentStiffness() const override { return true; }

    const char *giveInputRecordName() const override { return _IFT_IsoInterfaceDamageMaterial_Name; }
    const char *giveClassName() const override { return "IsoInterfaceDamageMaterial"; }

    FloatArrayF<3> giveEngTraction_3d(const FloatArrayF<3> &jump, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> give3dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    double computeEquivalentStrain(const FloatArrayF<3> &jump, GaussPoint *gp, TimeStep *tStep) const;

    /**
     * Computes the value of damage parameter omega, based on given value of equivalent strain.
     * @param kappa Equivalent strain measure.
     * @param strain Total strain vector in full form. (unnecessary?)
     * @param gp Integration point.
     * @return omega.
     */
    virtual double computeDamageParam(double kappa, const FloatArrayF<3> &strain, GaussPoint *gp) const;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<IsoInterfaceDamageMaterialStatus>(gp); }
};
} // end namespace oofem
#endif // isointerfacedamage01_h
