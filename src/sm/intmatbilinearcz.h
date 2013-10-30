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

/**
 *
 *	IntMatBilinearCZ
 *
 *	Bilinear cohesive zone model.
 *
 *  Created on: Oct 20, 2013
 *  @author: Erik Svenning
 */

#ifndef INTMATBILINEARCZ_H_
#define INTMATBILINEARCZ_H_

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

#include "dynamicinputrecord.h"

///@name Input fields for IntMatBilinearCZ
//@{
#define _IFT_IntMatBilinearCZ_Name "intmatbilinearcz"

#define _IFT_IntMatBilinearCZ_PenaltyStiffness "kn"
#define _IFT_IntMatBilinearCZ_g1c "g1c"
#define _IFT_IntMatBilinearCZ_g2c "g2c"
#define _IFT_IntMatBilinearCZ_mu "mu"
#define _IFT_IntMatBilinearCZ_gamma "gamma"
#define _IFT_IntMatBilinearCZ_sigf "sigf"
//@}

namespace oofem {

/**
 * This class implements associated Material Status for IntMatBilinearCZFagerstrom
 */
class IntMatBilinearCZStatus : public StructuralInterfaceMaterialStatus
{
public:
/*
    // material jump
    FloatArray oldMaterialJump;
    // temporary material jump
    FloatArray tempMaterialJump;
*/
    // damage variable
    double mDamageNew, mDamageOld;

    FloatArray mTractionOld, mTractionNew;
    FloatArray mJumpOld, mJumpNew;
/*
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
    FloatArray	alphav;
*/


public:

    IntMatBilinearCZStatus(int n, Domain *d, GaussPoint *g);

    virtual ~IntMatBilinearCZStatus();

//    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "IntMatBilinearCZStatus"; }
    virtual classType giveClassID() const { return MaterialStatusClass; }

//    double giveDamage() const { return mDamage; }
//    void setDamage(const double &iDamage) {mDamage = iDamage;}


/*
    double giveTempDamage() { return tempDamage; }

    const FloatArray &giveOldMaterialJump() {return oldMaterialJump; }
    const FloatArray &giveTempMaterialJump() {return tempMaterialJump; }

    const FloatArray &giveEffectiveMandelTraction() { return  QEffective; }
    const FloatArray &giveTempEffectiveMandelTraction() {return tempQEffective; }

    const FloatMatrix &giveTempInverseDefGrad() {return tempFInv; }
    const FloatMatrix &giveTempRotationMatrix() {return tempRot; }
    const FloatMatrix &giveTempIep() {return Iep; }
    const FloatArray &giveTempAlphav() {return alphav; }

*/
//    void letTempDamageBe(double v) { tempDamage = v; }
/*
    void letTempEffectiveMandelTractionBe(const FloatArray &v) { tempQEffective = v; }
    void letTempMaterialJumpBe(const FloatArray &v) { tempMaterialJump = v; }

    void letTempInverseDefGradBe(const FloatMatrix &v) { tempFInv = v; }
    void letTempRotationMatrix(const FloatMatrix &v) { tempRot = v; }
    void letTempIepBe(const FloatMatrix &v) { Iep = v; }
    void letTempAlphavBe(const FloatArray &v) { alphav = v; }

*/

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);


    //virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    //virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};



class IntMatBilinearCZ : public StructuralInterfaceMaterial {
public:
	IntMatBilinearCZ(int n, Domain *d);
	virtual ~IntMatBilinearCZ();

protected:
    /// Material parameters
    double mPenaltyStiffness;
    double mGIc;   // fracture energy, mode 1
    double mGIIc;  // fracture energy, mode 1
    double mSigmaF;  // max stress

    double mMu;    // loading function parameter
    double mGamma; // loading function parameter


    virtual int checkConsistency();
//    void give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
 //                                                                    GaussPoint *gp, TimeStep *atTime);

public:

    virtual int hasNonLinearBehaviour()   { return 1; }

//    virtual int hasMaterialModeCapability(MaterialMode mode); // remove
    virtual const char *giveClassName() const { return "IntMatBilinearCZ"; }
    virtual classType giveClassID() const { return MaterialStatusClass; }
    virtual const char *giveInputRecordName() const { return _IFT_IntMatBilinearCZ_Name; }


    virtual void giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                         const FloatMatrix &F, TimeStep *tStep);
    virtual void give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

private:
    double computeYieldFunction(const double &iTractionNormal, const double &iTractionTang) const;
    void computeTraction( FloatArray &oT, const FloatArray &iTTrial, const double &iPlastMultInc ) const;

public:
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new IntMatBilinearCZStatus(1, domain, gp); }
    virtual void printYourself();


};

} /* namespace oofem */
#endif /* INTMATBILINEARCZ_H_ */
