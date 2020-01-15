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

#ifndef rankinemat_h
#define rankinemat_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"

// this turns on or off a bunch of internal variables
// that allow tracking the distribution of dissipated energy
// (can be turned off if such information is not needed)
#define keep_track_of_dissipated_energy
//#undefine keep_track_of_dissipated_energy

///@name Input fields for RankineMat
//@{
#define _IFT_RankineMat_Name "rankmat"
#define _IFT_RankineMat_sig0 "sig0"
#define _IFT_RankineMat_h "h"
#define _IFT_RankineMat_a "a"
#define _IFT_RankineMat_plasthardtype "plasthardtype"
#define _IFT_RankineMat_delsigy "delsigy"
#define _IFT_RankineMat_yieldtol "yieldtol"
#define _IFT_RankineMat_gf "gf"
#define _IFT_RankineMat_ep "ep"
#define _IFT_RankineMat_damlaw "damlaw"
#define _IFT_RankineMat_param1 "param1"
#define _IFT_RankineMat_param2 "param2"
#define _IFT_RankineMat_param3 "param3"
#define _IFT_RankineMat_param4 "param4"
#define _IFT_RankineMat_param5 "param5"
//@}

namespace oofem {
class RankineMatStatus;

/**
 * This class implements an isotropic elastoplastic material
 * with Rankine yield condition, associated flow rule
 * and linear isotropic softening, and with isotropic damage
 * that leads to softening.
 *
 * It differs from other similar materials (such as RankinePlasticMaterial)
 * by implementation - here we use an efficient algorithm for this
 * specific model and implement the plane stress case and the vertex return.
 * Also, hardening and softening is incorporated.
 *
 */
class RankineMat : public StructuralMaterial
{
protected:
    /// Reference to the basic elastic material.
    LinearElasticMaterial *linearElasticMaterial = nullptr;

    /// Young's modulus.
    double E = 0.;

    /// Poisson's ratio.
    double nu = 0.;

    /// Initial hardening modulus.
    double H0 = 0.;

    /// Type of plastic hardening (0=linear, 1=exponential)
    int plasthardtype = 0;

    /// Final increment of yield stress (at infinite cumulative plastic strain)
    double delSigY = 0.;

    /// Initial (uniaxial) yield stress.
    double sig0 = 0.;

    /// Relative tolerance in yield condition
    double yieldtol = 0.;

    /// Parameter that controls damage evolution (a=0 turns damage off).
    double a = 0.;

    /// Total strain at peak stress sig0--Used only if plasthardtype=2
    double ep = 0.;

    /// Exponent in hardening law--Used only if plasthardtype=2
    double md = 0.;

    /// type of damage law (0=exponential, 1=exponential and  damage starts after peak stress sig0)
    int damlaw = 1;

    /// coefficient required when damlaw=1 or 2
    double param1 = 0.;

    /// coefficient required when  damlaw=1 or 2
    double param2 = 0.;

    /// coefficient required when damlaw=2
    double param3 = 0.;

    /// coefficient required when damlaw=2
    double param4 = 0.;

    /// coefficient required when damlaw=2
    double param5 = 0.;


public:
    RankineMat(int n, Domain * d);

    double evalYieldFunction(const FloatArray &sigPrinc, const double kappa) const;
    double evalYieldStress(const double kappa) const;
    double evalPlasticModulus(const double kappa) const;
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain) const;
    double computeDamage(GaussPoint *gp, TimeStep *tStep) const;
    double computeDamageParam(double tempKappa) const;
    double computeDamageParamPrime(double tempKappa) const;
    virtual double computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const;

    bool hasMaterialModeCapability(MaterialMode mode) const override;

    void initializeFrom(InputRecord &ir) override;

    // identification and auxiliary functions
    const char *giveInputRecordName() const override { return _IFT_RankineMat_Name; }
    const char *giveClassName() const override { return "RankineMat"; }

    /// Returns a reference to the basic elastic material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return ( a == 0. ); }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF<3> &reducesStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &totalStrain, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<1,1> give1dStressStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

protected:
    /**
     * Executive method used by local and gradient version.
     * (with different parameters gprime)
     */
    FloatMatrixF<3,3> evaluatePlaneStressStiffMtrx(MatResponseMode mode, GaussPoint *gp,
                                                   TimeStep *tStep, double gprime) const;

    /// Computes derivatives of final kappa with respect to final strain.
    void computeEta(FloatArray &answer, RankineMatStatus *status);

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
};

//=============================================================================


class RankineMatStatus : public StructuralMaterialStatus
{
protected:
    /// Plastic strain (initial).
    FloatArray plasticStrain;

    /// Plastic strain (final).
    FloatArray tempPlasticStrain;

    /// Effective stress (initial).
    FloatArray effStress;

    /// Effective stress (final).
    FloatArray tempEffStress;

    /// Cumulative plastic strain (initial).
    double kappa = 0.;

    /// Cumulative plastic strain (final).
    double tempKappa = 0.;

    /**
     * Increments of cumulative plastic strain
     * associated with the first and secomnd principal stress
     * (used in the case of vertex return, needed for stiffness)
     */
    double dKappa1 = 0., dKappa2 = 0.;

    /// Damage (initial).
    double damage = 0.;

    /// Damage (final).
    double tempDamage = 0.;

    /// Tangent shear stiffness (needed for tangent matrix).
    double tanG = 0.;

#ifdef keep_track_of_dissipated_energy
    /// Density of total work done by stresses on strain increments.
    double stressWork = 0.;
    /// Non-equilibrated density of total work done by stresses on strain increments.
    double tempStressWork = 0.;
    /// Density of dissipated work.
    double dissWork = 0.;
    /// Non-equilibrated density of dissipated work.
    double tempDissWork = 0.;
#endif

public:
    RankineMatStatus(GaussPoint * g);

    const FloatArray & givePlasticStrain() const { return plasticStrain; }

    double giveDamage() const { return damage; }
    double giveTempDamage() const { return tempDamage; }

    double giveCumulativePlasticStrain() const { return kappa; }
    double giveTempCumulativePlasticStrain() const { return tempKappa; }

    double giveDKappa(int i) const
    {
        if ( i == 1 ) {
            return dKappa1;
        } else {
            return dKappa2;
        }
    }

    double giveTangentShearStiffness()
    { return tanG; }

    const FloatArray &giveEffectiveStress() const { return effStress; }
    const FloatArray &giveTempEffectiveStress() const { return tempEffStress; }

    void letTempPlasticStrainBe(FloatArray values) { tempPlasticStrain = std :: move(values); }

    void letEffectiveStressBe(FloatArray values) { effStress = std :: move(values); }

    void letTempEffectiveStressBe(FloatArray values) { tempEffStress = std :: move(values); }

    void setTempCumulativePlasticStrain(double value) { tempKappa = value; }

    void setDKappa(double val1, double val2) {
        dKappa1 = val1;
        dKappa2 = val2;
    }

    void setTempDamage(double value) { tempDamage = value; }

    void setTangentShearStiffness(double value) { tanG = value; }

    const FloatArray &givePlasDef() { return plasticStrain; }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

#ifdef keep_track_of_dissipated_energy
    /// Returns the density of total work of stress on strain increments.
    double giveStressWork() { return stressWork; }
    /// Returns the temp density of total work of stress on strain increments.
    double giveTempStressWork() { return tempStressWork; }
    /// Sets the density of total work of stress on strain increments to given value.
    void setTempStressWork(double w) { tempStressWork = w; }
    /// Returns the density of dissipated work.
    double giveDissWork() { return dissWork; }
    /// Returns the density of temp dissipated work.
    double giveTempDissWork() { return tempDissWork; }
    /// Sets the density of dissipated work to given value.
    void setTempDissWork(double w) { tempDissWork = w; }
    /**
     * Computes the increment of total stress work and of dissipated work
     * (gf is the dissipation density per unit volume at complete failure,
     * it is needed only to determine which extremely small dissipation
     * can be set to zero to get clean results, but parameter gf can be
     * set to zero if not available).
     */
    void computeWork_PlaneStress(GaussPoint *gp, double gf);
    void computeWork_1d(GaussPoint *gp, double gf);
#endif

    const char *giveClassName() const override { return "RankineMatStatus"; }
};
} // end namespace oofem
#endif // rankinemat_h
