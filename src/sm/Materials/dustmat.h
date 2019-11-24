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

#ifndef dustmat_h
#define dustmat_h

#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/isolinearelasticmaterial.h"

///@name Input fields for DustMaterial
//@{
#define _IFT_DustMaterial_Name "dustmat"
#define _IFT_DustMaterial_alpha "alpha"
#define _IFT_DustMaterial_beta "beta"
#define _IFT_DustMaterial_lambda "lambda"
#define _IFT_DustMaterial_theta "theta"
#define _IFT_DustMaterial_ft "ft"
#define _IFT_DustMaterial_hardeningType "hardeningtype"
#define _IFT_DustMaterial_mStiff "mstiff"
#define _IFT_DustMaterial_rEll "rell"
#define _IFT_DustMaterial_x0 "x0"
#define _IFT_DustMaterial_newtonTol "newtontol"
#define _IFT_DustMaterial_newtonIter "newtoniter"
#define _IFT_DustMaterial_wHard "whard"
#define _IFT_DustMaterial_dHard "dhard"
//@}

namespace oofem {
/**
 * This class implements material status for dust material model
 * @see DustMaterial
 * @author Jan Stransky
 */
class DustMaterialStatus : public StructuralMaterialStatus
{
public:
    /// Values of history variable stateFlag
    enum stateFlagValues {
        DM_Elastic,
        DM_Unloading,
        DM_Yielding1,
        DM_Yielding2,
        DM_Yielding3
    };

protected:
    /// Plastic strain
    FloatArrayF<6> plasticStrain;
    FloatArrayF<6> tempPlasticStrain;

    /// Hardening parameter q
    double q = 0;
    double tempQ = 0;

    /// Current bulk modulus
    double bulkModulus = 0;
    /// Current shear modulus
    double shearModulus = 0;
    /// Current Young's modulus
    double youngsModulus = 0;

    /// Indicates the state (i.e. elastic, yielding, unloading) of the Gauss point
    int stateFlag = 0;
    int tempStateFlag = 0;

public:
    /// Constructor
    DustMaterialStatus(GaussPoint * gp, double q0);

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const char *giveClassName() const override { return "DustMaterialStatus"; }

    /**
     * Get the full plastic strain vector from the material status.
     * @return Plastic strain.
     */
    const FloatArrayF<6> &givePlasticStrain() const { return plasticStrain; }
    /**
     * Get the full plastic strain vector from the material status.
     * @return Volumetric part of plastic strain.
     */
    double giveVolumetricPlasticStrain() const { return ( plasticStrain.at(1) + plasticStrain.at(2) + plasticStrain.at(3) ) / 3.0; }
    /**
     * Get the value of hardening variable q from the material status.
     * @return Hardenenig variable q
     */
    double giveQ() const { return q; }
    /**
     * Get the state flag from the material status.
     * @return State flag (i.e. elastic, unloading, yielding, vertex case yielding)
     */
    int giveStateFlag() const { return stateFlag; }

    /**
     * Get the temp value of the full plastic strain vector from the material status.
     * @return Temp value of plastic strain vector.
     */
    const FloatArrayF<6> &giveTempPlasticStrain() const { return tempPlasticStrain; }
    /**
     * Get the temp value of the hardening variable q from the material status.
     * @return Temp value of hardening variable kappaP.
     */
    double giveTempQ() const { return tempQ; }
    /**
     * Get the temp value of the state flag from the material status.
     * @return The temp value of the state flag.
     */
    int giveTempStateFlag() const { return tempStateFlag; }

    /**
     * Assign the temp value of plastic strain.
     * @param v New temp value o plastic strain.
     */
    void letTempPlasticStrainBe(const FloatArray &v) { tempPlasticStrain = v; }
    /**
     * Assign the temp value of variable q
     * @param v New temp value of variable q
     */
    void letTempQBe(double v) { tempQ = v; }
    /**
     * Assign the temp value of the state flag.
     * @param v New temp value of the state flag.
     */
    void letTempStateFlagBe(int v) { tempStateFlag = v; }

    /**
     * Assign the value of plastic strain.
     * @param v New value of plastic strain.
     */
    void letPlasticStrainBe(const FloatArrayF<6> &v) { plasticStrain = v; }
    /**
     * Assign the value of variable q
     * @param v New value of variable q
     */
    void letQBe(double v) { q = v; }

    /**
     * Assign the value of actual bulk modulus of the status
     * @param v New value of bulk modulus
     */
    void setBulkModulus(double v) { bulkModulus = v; }
    /**
     * Assign the value of actual shear modulus of the status
     * @param v New value of shear modulus
     */
    void setShearModulus(double v) { shearModulus = v; }
    /**
     * Assign the value of actual Young's modulus of the status
     * @param v New value of Young's modulus
     */
    void setYoungsModulus(double v) { youngsModulus = v; }
    /**
     * Get the value of actual bulk modulus of the status
     * @return Value of bulk modulus
     */
    double giveBulkModulus() const { return bulkModulus; }
    /**
     * Get the value of actual shear modulus of the status
     * @return Value of shear modulus
     */
    double giveShearModulus() const { return shearModulus; }
    /**
     * Get the value of actual Young's modulus of the status
     * @return Value of Young's modulus
     */
    double giveYoungsModulus() const { return youngsModulus; }
};

/**
 * This class implements nonassociated multisurface plasticity model.
 * Yield function depends on I1 and J2 invariants of stress tensor and
 * internal variable q.
 * The plasticity surface consists of three surfaces.
 * Shear dominant surface (f1) is non linear Drucker Prager cone.
 * Compressive dominant part (f2) is elliptic cap over the cone.
 * Tension dominant part (f3) is plane, cutting the cone perpendicularly to
 * hydrostatic axis.
 * f1 and f3 are constant. The position of f2 surface is determined by internal
 * variable q. During plastic yielding f1 or f3 (shear or tension dominant state),
 * q is increasing and the admissible elastic domain is shrinking. During
 * f2 plastic yielding (compressive dominant), the parameter q is decreasing and the
 * admissible elastic domain is becoming larger.
 *
 * This model is described in PhD thesis of Tobias Erhard
 * Strategien zur numerischen Modellierung transierten Impaktvorgange bei nichtlinearem Materialverhalten
 * Universitat Stuttgart 2004
 *
 * @author Jan Stransky
 */
class DustMaterial : public StructuralMaterial
{
protected:
    /// Pointer for linear elastic material
    IsotropicLinearElasticMaterial LEMaterial;

    /// Parameter determining shape of yield surface
    double alpha = 0.;
    /// Parameter determining shape of yield surface
    double beta = 0.;
    /// Parameter determining shape of yield surface
    double lambda = 0.;
    /// Parameter determining shape of yield surface
    double theta = 0.;
    /// Parameter determining shape of yield surface (param T in original publication)
    double ft = 0.;
    /// Parameter determining shape of yield surface (param R in original publication)
    double rEll = 0.;
    /// Parameter determining hardening type
    int hardeningType = 0;
    /// Parameter increasing stiffness (parameter M in original publication)
    double mStiff = 0.;
    /// Parameter determining shape of yield surface (param X0 in original publication)
    double x0 = 0.;
    /// Parameter determining shape of yield surface
    double q0 = 0.;
    /// Parameter determining hardening law (parameter W in original publication)
    double wHard = 0.;
    /// Parameter determining hardening law (parameter D in original publication)
    double dHard = 0.;
    /// Tollerance for iterative methods
    double newtonTol = 0.;
    /// Maximum number of iterations for iterative methods
    int newtonIter = 0;

    /**
     * Auxiliary equation Fe (7.8)
     * @param i1 Trace of stress tensor
     * @return Fe
     */
    double functionFe(double i1) const;
    /**
     * Derivative by i1 of auxiliary equation (7.8)
     * @param i1 Trace of stress tensor
     * @return @f$ \frac{\partial Fe}{\partial i_1 } @f$
     */
    double functionFeDI1(double i1) const;
    /**
     * Second derivative by i1 of auxiliary equation (7.8)
     * @param i1 Trace of stress tensor
     * @return @f$ \frac{\partial^2 Fe}{\partial i1^2 } @f$
     */
    double functionFeDI1DI1(double i1) const;
    /**
     * Auxiliary equation Fc (7.9)
     * @param i1 Trace of stress tensor
     * @param rho Second Haigh-Westergaard coordinate
     * @param q Parameter q
     * @return Fc
     */
    double functionFc(double rho, double i1, double q) const;
    /**
     * Auxiliary equation X (7.11)
     * @param q Parameter q
     * @return X
     */
    double functionX(double q) const;
    /**
     * Derivative by q of auxiliary equation X (7.11)
     * @param q Parameter q
     * @return @f$ \frac{\partial X }{\partial q} @f$
     */
    double functionXDQ(double q) const;
    /**
     * Yield function 1 (shear dominant), equation 7.5
     * @param rho Second Haigh-Westergaard coordinate
     * @param i1 Trace of stress tensor
     * @return Value
     */
    double yieldFunction1(double rho, double i1) const;
    /**
     * Yield function 2 (compression dominant), equation 7.6
     * @param rho Second Haigh-Westergaard coordinate
     * @param i1 Trace of stress tensor
     * @param q Parameter q
     * @return value
     */
    double yieldFunction2(double rho, double i1, double q) const;
    /**
     * Yield function 3 (tension dominant), equation 7.7
     * @param i1 Trace of stress tensor
     * @return value
     */
    double yieldFunction3(double i1) const;
    /**
     * Solves q0 according to given parameters, equation 7.12
     * @param answer Result
     */
    void solveQ0(double &answer) const;
    /**
     * Computes and sets all elastic moduli, with possible stiffening. Equation 7.4
     * @param bulkModulus Bulk modulus
     * @param shearModulus Shear modulus
     * @param gp Gauss point
     */
    void computeAndSetBulkAndShearModuli(double &bulkModulus, double &shearModulus, GaussPoint *gp) const;
    /**
     * Perform stress return and update all internal variables
     * @param gp Gauss point
     * @param strain Strain
     */
    void performStressReturn(GaussPoint *gp, const FloatArrayF<6> &strain) const;
    /**
     * Computes direction of plastic yielding m1, equation 7.17
     * @param stressDeviator Deviator of stress tensor
     * @param rho Second Haigh-Westergaard coordinate
     * @param i1 Trace of stress tensor
     * @param q Parameter q
     * @param return Result
     */
    FloatArrayF<6> computePlastStrainDirM1(const FloatArrayF<6> &stressDeviator, double rho, double i1, double q) const;
    /**
     * Computes direction of plastic yielding m2, equation 7.19
     * @param stressDeviator Deviator of stress tensor
     * @param rho Second Haigh-Westergaard coordinate
     * @param i1 Trace of stress tensor
     * @param q Parameter q
     * @param return Result
     */
    FloatArrayF<6> computePlastStrainDirM2(const FloatArrayF<6> &stressDeviator, double rho, double i1, double q) const;
    /**
     * Computes direction of plastic yielding m2, equation 7.18
     * @param stressDeviator Deviator of stress tensor
     * @param rho Second Haigh-Westergaard coordinate
     * @param i1 Trace of stress tensor
     * @param q Parameter q
     * @param return Result
     */
    FloatArrayF<6> computePlastStrainDirM3(const FloatArrayF<6> &stressDeviator, double rho, double i1, double q) const;
    /**
     * Auxiliary equation H (7.33 or 7.34)
     * @param q Parameter q from previous step
     * @param tempQ Parameter tempQ
     * @return H
     */
    double functionH(double q, double tempQ) const;
    /**
     * Derivative by tempQ of auxiliary equation H (7.33 or 7.34)
     * @param tempQ Parameter tempQ
     * @return @f$ \frac{\partial X}{\partial tempQ} @f$
     */
    double functionHDQ(double tempQ) const;
    /**
     * Auxiliary equation I1 (7.32)
     * @param q Parameter q from previous step
     * @param tempQ Parameter tempQ
     * @param i1 Trace of stress tensor
     * @param bulkModulus Bulk modulus
     * @return I1
     */
    double functionI1(double q, double tempQ, double i1, double bulkModulus) const;
    /**
     * Derivative by tempQ of auxiliary equation I1 (7.32)
     * @param tempQ Parameter tempQ
     * @param bulkModulus Bulk modulus
     * @return @f$ \frac{\partial I1}{\partial tempQ} @f$
     */
    double functionI1DQ(double tempQ, double bulkModulus) const;
    /**
     * Performs stress return of case of yield function F1, computes new value of tempQ and sets it to status. Equation 7.31
     * @param i1 Trace of stress tensor
     * @param rho Second Haigh-Westergaard coordinate
     * @param gp Gauss point
     */
    void performF1return(double i1, double rho, GaussPoint *gp) const;
    /**
     * Performs stress return of case of yield function F2, computes new value of tempQ and sets it to status. Equation 7.38
     * @param i1 Trace of stress tensor
     * @param rho Second Haigh-Westergaard coordinate
     * @param gp Gauss point
     */
    void performF2return(double i1, double rho, GaussPoint *gp) const;
    /**
     * Computes tempQ from volumetric plastic strain increment, equation 7.44
     * @param answer Result tempQ
     * @param q Parameter q from previous step
     * @param deltaVolumetricPlasticStrain Volumetric plastic strain increment
     */
    void computeQFromPlastVolEps(double &answer, double q, double deltaVolumetricPlasticStrain) const;
    /**
     * Computed value of plastic multiplier for F2 yield function, equation 7.39
     * @param tempQ Parameter tempQ
     * @param q Parameter q from previous step
     * @param i1 Trace of stress tensor
     * @param bulkModulus Bulk modulus
     */
    double computeDeltaGamma2(double tempQ, double q, double i1, double bulkModulus) const;
    /**
     * Computed derivative by tempQ of equation 7.39
     * @param tempQ Parameter tempQ
     * @param q Parameter q from previous step
     * @param i1 Trace of stress tensor
     * @param bulkModulus Bulk modulus
     */
    double computeDeltaGamma2DQ(double tempQ, double q, double i1, double bulkModulus) const;
    /**
     * Equation 7.38
     * @param tempQ Parameter tempQ
     * @param q Parameter q from previous step
     * @param i1 Trace of stress tensor
     * @param rho Second Haigh-Westergaard coordinate
     * @param bulkModulus Bulk modulus
     * @param shearModulus Shear modulus
     */
    double fTempR2(double tempQ, double q, double i1, double rho, double bulkModulus, double shearModulus) const;

public:
    /// Constructor
    DustMaterial(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;

    const char *giveClassName() const override { return "DustMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_DustMaterial_Name; }

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;

    int setIPValue(const FloatArray &value, GaussPoint *gp, InternalStateType type) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    FloatArrayF<6> giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override
    {
        return LEMaterial.giveThermalDilatationVector(gp, tStep);
    }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    double giveQ0() const { return q0; }
};
} // end namespace oofem
#endif // dustmat_h
