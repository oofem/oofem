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
//@}

namespace oofem {
/**
 * This class implements associated Material Status for IntMatBilinearCZJansson
 */
class IntMatBilinearCZJanssonStatus : public StructuralInterfaceMaterialStatus
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
    IntMatBilinearCZJanssonStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~IntMatBilinearCZJanssonStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "IntMatBilinearCZJanssonStatus"; }

    double giveDamage() { return damage; }
    double giveTempDamage() { return tempDamage; }
    bool giveOldDamageDev() {return oldDamageDev;}

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
    void letTempEffectiveMandelTractionBe(FloatArray v) { tempQEffective = std :: move(v); }
    void letTempMaterialJumpBe(FloatArray v) { tempMaterialJump = std :: move(v); }
    void letTempDamageDevBe(bool v) { tempDamageDev = v; }
    void letOldDamageDevBe(bool v) { oldDamageDev = v; }

    void letTempdTdJBe(FloatMatrix &v) { temp_dTdJ = v; }

    void letTempInverseDefGradBe(FloatMatrix v) { tempFInv = std :: move(v); }
    void letTempRotationMatrix(FloatMatrix v) { tempRot = std :: move(v); }
    void letTempIepBe(FloatMatrix v) { Iep = std :: move(v); }
    void letTempAlphavBe(FloatArray v) { alphav = std :: move(v); }



    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    //virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    //virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
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
    /// Material parameters
    double kn0;   // initial normal stiffness
    double ks0;   // initial shear stiffness
    double knc;   // stiffness in compression
    double GIc;   // fracture energy, mode 1
    double GIIc;  // fracture energy, mode 1
    double sigf;  // max stress

    double mu;    // loading function parameter
    double gamma; // loading function parameter


    virtual int checkConsistency();
    void give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                                GaussPoint *gp, TimeStep *tStep);

public:
    /// Constructor
    IntMatBilinearCZJansson(int n, Domain * d);
    /// Destructor
    virtual ~IntMatBilinearCZJansson();

    virtual int hasNonLinearBehaviour()   { return 1; }

    virtual const char *giveClassName() const { return "IntMatBilinearCZJansson"; }
    virtual const char *giveInputRecordName() const { return _IFT_IntMatBilinearCZJansson_Name; }


    virtual void giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                        const FloatMatrix &F, TimeStep *tStep);
    virtual void give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual bool hasAnalyticalTangentStiffness() const { return true; }

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new IntMatBilinearCZJanssonStatus(1, domain, gp); } //@Martin: Why new?
    void printYourself();
protected:
};
} // end namespace oofem
#endif // isointerfacedamage01_h
