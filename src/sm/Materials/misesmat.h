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

#ifndef misesmat_h
#define misesmat_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "scalarfunction.h"

///@name Input fields for MisesMat
//@{
#define _IFT_MisesMat_Name "misesmat"
#define _IFT_MisesMat_sig0 "sig0"
#define _IFT_MisesMat_h "h"
#define _IFT_MisesMat_omega_crit "omega_crit"
#define _IFT_MisesMat_a "a"
#define _IFT_MisesMat_yieldTol "yieldtol"
//@}

namespace oofem {
class GaussPoint;
class Domain;

/**
 * This class implements an isotropic elastoplastic material
 * with Mises yield condition, associated flow rule
 * and linear isotropic hardening.
 *
 * It differs from other similar materials (such as J2Mat)
 * by implementation - here we use the radial return, which
 * is the most efficient algorithm for this specific model.
 * Also, an extension to large strain will be available.
 *
 * The model also exemplifies how to implement non-3d material modes, in this case 1D,
 * by overloading the default implementations that iterates over the 3D method.
 */
class MisesMat : public StructuralMaterial
{
protected:
    /// Reference to the basic elastic material.
    IsotropicLinearElasticMaterial linearElasticMaterial;

    /// Elastic shear modulus.
    double G;

    /// Elastic bulk modulus.
    double K;

    /// Hardening modulus.
    double H;

    /// Initial (uniaxial) yield stress.
    ScalarFunction sig0;

    /// critical(maximal) damage.
    double omega_crit;
    /// exponent in damage function.
    double a;

    /// tolerance for the yield function in RRM algorithm.
    double yieldTol;

public:
    MisesMat(int n, Domain *d);
    virtual ~MisesMat() {}

    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep);
    void performPlasticityReturn_PlaneStress(GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep);

    double computeYieldStress(double kappa, GaussPoint *gp, TimeStep *tStep);
    double computeYieldStressPrime(double kappa);

    double computeDamage(GaussPoint *gp, TimeStep *tStep);
    double computeDamageParam(double tempKappa);
    double computeDamageParamPrime(double tempKappa);
    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);

    IRResultType initializeFrom(InputRecord *ir) override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    const char *giveInputRecordName() const override { return _IFT_MisesMat_Name; }
    const char *giveClassName() const override { return "MisesMat"; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                       MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *tStep) override;


    void givePlaneStressStiffMtrx(FloatMatrix &answer,
                                  MatResponseMode mmode, GaussPoint *gp,
                                  TimeStep *tStep) override;


    void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep) override;

    void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep) override;

    void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep) override;


    double give(int aProperty, GaussPoint *gp, TimeStep *tStep);
    double giveTemperature(GaussPoint *gp, TimeStep *tStep);

protected:
    void computeGLPlasticStrain(const FloatMatrix &F, FloatMatrix &Ep, FloatMatrix b, double J);


    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
};

//=============================================================================


class MisesMatStatus : public StructuralMaterialStatus
{
protected:
    /// Plastic strain (initial).
    FloatArray plasticStrain;

    /// Plastic strain (final).
    FloatArray tempPlasticStrain;

    /// Deviatoric trial stress - needed for tangent stiffness.
    FloatArray trialStressD;

    /// volumetric trial stress - needed for tangent stiffness.
    double trialStressV;

    FloatArray effStress;
    FloatArray tempEffStress;

    /// Cumulative plastic strain (initial).
    double kappa;

    /// Cumulative plastic strain (final).
    double tempKappa;

    /// damage variable (initial).
    double damage;

    /// damage variable (final).
    double tempDamage;


public:
    MisesMatStatus(GaussPoint *g);

    const FloatArray &givePlasticStrain() const { return plasticStrain; }

    const FloatArray &giveTrialStressDev() const { return trialStressD; }

    double giveTrialStressVol() const { return trialStressV; }

    double giveDamage() const { return damage; }
    double giveTempDamage() const { return tempDamage; }

    double giveCumulativePlasticStrain() const { return kappa; }
    double giveTempCumulativePlasticStrain() const { return tempKappa; }

    const FloatArray &giveTempEffectiveStress() const { return tempEffStress; }
    const FloatArray &giveEffectiveStress() const { return effStress; }

    void letTempPlasticStrainBe(const FloatArray &values) { tempPlasticStrain = values; }
    const FloatArray &getTempPlasticStrain() const { return tempPlasticStrain; }

    void letTrialStressDevBe(const FloatArray &values) { trialStressD = values; }

    void letEffectiveStressBe(const FloatArray &values) { effStress = values; }

    void letTempEffectiveStressBe(const FloatArray &values) { tempEffStress = values; }


    void setTrialStressVol(double value) { trialStressV = value; }




    void setTempCumulativePlasticStrain(double value) { tempKappa = value; }

    void setTempDamage(double value) { tempDamage = value; }

    const FloatArray &givePlasDef() { return plasticStrain; }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;

    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const char *giveClassName() const override { return "MisesMatStatus"; }
};
} // end namespace oofem
#endif // misesmat_h
