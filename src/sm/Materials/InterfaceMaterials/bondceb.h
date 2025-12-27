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

#ifndef bondceb_h
#define bondceb_h

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for BondCEBMaterial
//@{
#define _IFT_BondCEBMaterial_Name "bondceb"
#define _IFT_BondCEBMaterial_kn "kn"
#define _IFT_BondCEBMaterial_ks "ks"
#define _IFT_BondCEBMaterial_s1 "s1"
#define _IFT_BondCEBMaterial_s2 "s2"
#define _IFT_BondCEBMaterial_s3 "s3"
#define _IFT_BondCEBMaterial_al "al"
#define _IFT_BondCEBMaterial_taumax "taumax"
#define _IFT_BondCEBMaterial_tauf "tauf"
//@}

namespace oofem {

/**
 * This class implements associated status to BondCEBInterfaceMaterial.
 */
class BondCEBMaterialStatus : public StructuralInterfaceMaterialStatus
{
protected:
    /// Cumulative slip.
    double kappa = 0.;
    /// Non-equilibrated cumulative slip.
    double tempKappa = 0.;

public:
    /// Constructor
    BondCEBMaterialStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    /// Returns the last equilibrated cumulative slip.
    double giveKappa() const { return kappa; }
    /// Returns the temporary cumulative slip.
    double giveTempKappa() const { return tempKappa; }
    /// Sets the temporary cumulative slip to the given value.
    void setTempKappa(double newKappa) { tempKappa = newKappa; }

    const char *giveClassName() const override { return "BondCEBMaterialStatus"; }

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
class BondCEBMaterial : public StructuralInterfaceMaterial
{
protected:
    /// Normal elastic stiffness.
    double kn = 0.;
    /// Shear elastic stiffness.
    double ks = 0.;
    /// Shear strength.
    double taumax = 0.;
    /// Residual shear stress.
    double tauf = 0.;
    /// Characteristic slip values.
    double s0 = 0., s1 = 0., s2 = 0., s3 = 0.;
    /// Exponent.
    double alpha = 0.4;

public:
    /// Constructor
    BondCEBMaterial(int n, Domain * d);

    bool hasAnalyticalTangentStiffness() const override { return true; }

    const char *giveInputRecordName() const override { return _IFT_BondCEBMaterial_Name; }
    const char *giveClassName() const override { return "BondCEBMaterial"; }

    FloatArrayF<3> giveEngTraction_3d(const FloatArrayF<3> &jump, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> give3dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<BondCEBMaterialStatus>(gp); }

protected:
    double evaluateBondStress(const double kappa) const;
};
} // end namespace oofem
#endif // bondceb_h
