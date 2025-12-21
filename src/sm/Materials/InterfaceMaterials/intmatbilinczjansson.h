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

#ifndef intmatbilinearczmaterialJansson_h
#define intmatbilinearczmaterialJansson_h

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for IntMatBilinearCZJansson
//@{
#define _IFT_IntMatBilinearCZJansson_Name "intmatbilinearczjansson"
#define _IFT_IntMatBilinearCZJansson_kn "kn"
#define _IFT_IntMatBilinearCZJansson_ks "ks"
#define _IFT_IntMatBilinearCZJansson_knc "knc"
#define _IFT_IntMatBilinearCZJansson_g1c "g1c"
#define _IFT_IntMatBilinearCZJansson_g2c "g2c"
#define _IFT_IntMatBilinearCZJansson_mu "mu"
#define _IFT_IntMatBilinearCZJansson_gamma "gamma"
#define _IFT_IntMatBilinearCZJansson_sigf "sigf"
#define _IFT_IntMatBilinearCZJansson_semiexplicit "semiexplicit"
//@}

namespace oofem {
/**
 * This class implements associated Material Status for IntMatBilinearCZJansson
 */
class IntMatBilinearCZJanssonStatus : public StructuralInterfaceMaterialStatus
{
protected:
    // material jump
    FloatArrayF<3> oldMaterialJump;
    // temporary material jump
    FloatArrayF<3> tempMaterialJump;

    // damage variable
    double damage = 0.;
    // temporary damage value
    double tempDamage = 0.;

    // Effective Mandel traction
    FloatArrayF<3> QEffective;
    // Temporary effective Mandel traction
    FloatArrayF<3> tempQEffective;

    // Temporary inverse of deformation gradient
    FloatMatrixF<3,3> tempFInv;

    // Temporary array for coordinate transformation
    FloatMatrixF<3,3> tempRot;

    // tempArrays for stiffness calculation
    FloatMatrixF<3,3> Iep;
    FloatArrayF<3> alphav;

    // indicator for davae development of preceding time step
    bool tempDamageDev = false;
    bool oldDamageDev = false;

    // tangent stiffness from previous time step
    FloatMatrixF<3,3> temp_dTdJ;
    FloatMatrixF<3,3> old_dTdJ;

public:
    /// Constructor
    IntMatBilinearCZJanssonStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "IntMatBilinearCZJanssonStatus"; }

    double giveDamage() const override { return damage; }
    double giveTempDamage() const override { return tempDamage; }
    bool giveOldDamageDev() const { return oldDamageDev; }

    const FloatArrayF<3> &giveOldMaterialJump() const { return oldMaterialJump; }
    const FloatArrayF<3> &giveTempMaterialJump() const { return tempMaterialJump; }

    const FloatArrayF<3> &giveEffectiveMandelTraction() const { return QEffective; }
    const FloatArrayF<3> &giveTempEffectiveMandelTraction() const { return tempQEffective; }

    const FloatMatrixF<3,3> &giveTempInverseDefGrad() const { return tempFInv; }
    const FloatMatrixF<3,3> &giveTempRotationMatrix() const { return tempRot; }
    const FloatMatrixF<3,3> &giveTempIep() const { return Iep; }
    const FloatArrayF<3> &giveTempAlphav() const { return alphav; }
    const FloatMatrixF<3,3> &giveOlddTdJ() const { return old_dTdJ; }


    void letTempDamageBe(double v) { tempDamage = v; }
    void letTempEffectiveMandelTractionBe(const FloatArrayF<3> &v) { tempQEffective = std :: move(v); }
    void letTempMaterialJumpBe(const FloatArrayF<3> &v) { tempMaterialJump = std :: move(v); }
    void letTempDamageDevBe(bool v) { tempDamageDev = v; }
    void letOldDamageDevBe(bool v) { oldDamageDev = v; }

    void letTempdTdJBe(const FloatMatrix &v) { temp_dTdJ = v; }

    void letTempInverseDefGradBe(const FloatMatrixF<3,3> &v) { tempFInv = std :: move(v); }
    void letTempRotationMatrix(const FloatMatrixF<3,3> &v) { tempRot = std :: move(v); }
    void letTempIepBe(const FloatMatrixF<3,3> &v) { Iep = std :: move(v); }
    void letTempAlphavBe(const FloatArrayF<3> &v) { alphav = std :: move(v); }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    //void saveContext(DataStream &stream, ContextMode mode) override;
    //void restoreContext(DataStream &stream, ContextMode mode) override;
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
class IntMatBilinearCZJansson : public StructuralInterfaceMaterial
{
protected:
    double kn0 = 0.;   // initial normal stiffness
    double ks0 = 0.;   // initial shear stiffness
    double knc = 0.;   // stiffness in compression
    double GIc = 0.;   // fracture energy, mode 1
    double GIIc = 0.;  // fracture energy, mode 1
    double sigf = 0.;  // max stress

    double mu = 0.;    // loading function parameter
    double gamma = 0.; // loading function parameter

    bool mSemiExplicit = false; // If semi-explicit time integration should be used

    int checkConsistency() override;
    void give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                                GaussPoint *gp, TimeStep *tStep);

public:
    IntMatBilinearCZJansson(int n, Domain * d);

    const char *giveClassName() const override { return "IntMatBilinearCZJansson"; }
    const char *giveInputRecordName() const override { return _IFT_IntMatBilinearCZJansson_Name; }

    FloatArrayF<3> giveFirstPKTraction_3d(const FloatArrayF<3> &jump, const FloatMatrixF<3,3> &F, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> give3dStiffnessMatrix_dTdj(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    bool hasAnalyticalTangentStiffness() const override { return true; }

    void initializeFrom(InputRecord &ir) override;

    FloatArray giveInterfaceStrength() override { return {this->sigf*this->gamma,this->sigf*this->gamma,this->sigf}; }

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<IntMatBilinearCZJanssonStatus>(gp); }
    void printYourself() override;
};
} // end namespace oofem
#endif // isointerfacedamage01_h
