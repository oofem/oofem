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

#ifndef intmatbilinearczmaterialFagerstrom_h
#define intmatbilinearczmaterialFagerstrom_h

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for IntMatBilinearCZFagerstrom
//@{
#define _IFT_IntMatBilinearCZFagerstrom_Name "intmatbilinearczfagerstrom"
#define _IFT_IntMatBilinearCZFagerstrom_kn "kn"
#define _IFT_IntMatBilinearCZFagerstrom_ks "ks"
#define _IFT_IntMatBilinearCZFagerstrom_knc "knc"
#define _IFT_IntMatBilinearCZFagerstrom_g1c "g1c"
#define _IFT_IntMatBilinearCZFagerstrom_g2c "g2c"
#define _IFT_IntMatBilinearCZFagerstrom_mu "mu"
#define _IFT_IntMatBilinearCZFagerstrom_gamma "gamma"
#define _IFT_IntMatBilinearCZFagerstrom_sigf "sigf"
//@}

namespace oofem {
/**
 * This class implements associated Material Status for IntMatBilinearCZFagerstrom
 */
class IntMatBilinearCZFagerstromStatus : public StructuralInterfaceMaterialStatus
{
protected:
    // material jump
    FloatArray oldMaterialJump;
    // temporary material jump
    FloatArray tempMaterialJump;

    // damage variable
    double damage;
    // temporary damage value
    double tempDamage;

    // Effective Mandel traction
    FloatArray QEffective;
    // Temporary effective Mandel traction
    FloatArray tempQEffective;

    // Temporary inverse of deformation gradient
    FloatMatrix tempFInv;

    // Temporary array for coordinate transformation
    FloatMatrix tempRot;

    // tempArrays for stiffness calculation
    FloatMatrix Iep;
    FloatArray alphav;

    // indicator for davae development of preceding time step
    bool tempDamageDev;
    bool oldDamageDev;

    // tangent stiffness from previous time step
    FloatMatrix temp_dTdJ;
    FloatMatrix old_dTdJ;

public:
    /// Constructor
    IntMatBilinearCZFagerstromStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~IntMatBilinearCZFagerstromStatus();

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    // definition
    const char *giveClassName() const override { return "IntMatBilinearCZFagerstromStatus"; }

    double giveDamage() override { return damage; }
    double giveTempDamage() override { return tempDamage; }
    bool giveOldDamageDev() { return oldDamageDev; }

    const FloatArray &giveOldMaterialJump() { return oldMaterialJump; }
    const FloatArray &giveTempMaterialJump() { return tempMaterialJump; }

    const FloatArray &giveEffectiveMandelTraction() { return QEffective; }
    const FloatArray &giveTempEffectiveMandelTraction() { return tempQEffective; }

    const FloatMatrix &giveTempInverseDefGrad() { return tempFInv; }
    const FloatMatrix &giveTempRotationMatrix() { return tempRot; }
    const FloatMatrix &giveTempIep() { return Iep; }
    const FloatArray &giveTempAlphav() { return alphav; }
    const FloatMatrix &giveOlddTdJ() {return old_dTdJ; }


    void letTempDamageBe(double v) { tempDamage = v; }
    void letTempDamageDevBe(bool v) { tempDamageDev = v; }
    void letOldDamageDevBe(bool v) { oldDamageDev = v; }
    void letTempEffectiveMandelTractionBe(FloatArray v) { tempQEffective = std :: move(v); }
    void letTempMaterialJumpBe(FloatArray v) { tempMaterialJump = std :: move(v); }

    void letTempdTdJBe(FloatMatrix &v) { temp_dTdJ = v; }

    void letTempInverseDefGradBe(FloatMatrix v) { tempFInv = std :: move(v); }
    void letTempRotationMatrix(FloatMatrix v) { tempRot = std :: move(v); }
    void letTempIepBe(FloatMatrix v) { Iep = std :: move(v); }
    void letTempAlphavBe(FloatArray v) { alphav = std :: move(v); }

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
class IntMatBilinearCZFagerstrom : public StructuralInterfaceMaterial
{
protected:
    /// Material parameters
    double kn0;   // initial normal stiffness
    double ks0;   // initial shear stiffness
    double knc;   // stiffness in compression
    double GIc;   // fracture energy, mode 1
    double GIIc;  // fracture energy, mode 1
    double sigf;  // max stress

    double mu;    // loading function parameter
    double gamma; // loading function parameter


    int checkConsistency() override;
    void give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                                GaussPoint *gp, TimeStep *tStep);

public:
    /// Constructor
    IntMatBilinearCZFagerstrom(int n, Domain * d);
    /// Destructor
    virtual ~IntMatBilinearCZFagerstrom();

    const char *giveClassName() const override { return "IntMatBilinearCZFagerstrom"; }
    const char *giveInputRecordName() const override { return _IFT_IntMatBilinearCZFagerstrom_Name; }

    void giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                const FloatMatrix &F, TimeStep *tStep) override;
    void give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

    /**
     * Tells if the model has implemented analytical tangent stiffness.
     * If not, the tangent must be computed numerically.
     */
    bool hasAnalyticalTangentStiffness() const override { return true; }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    FloatArray giveInterfaceStrength() override { return {this->sigf*this->gamma,this->sigf*this->gamma,this->sigf}; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new IntMatBilinearCZFagerstromStatus(1, domain, gp); } //@Martin: Why new?
    void printYourself() override;
};
} // end namespace oofem
#endif // isointerfacedamage01_h
