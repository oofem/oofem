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

#ifndef druckerpragercatmat_h
#define druckerpragercatmat_h

#include "Materials/ConcreteMaterials/mplasticmaterial2.h"

///@name Input fields for DruckerPragerCutMat
//@{
#define _IFT_DruckerPragerCutMat_Name "druckerpragercutmat"
#define _IFT_DruckerPragerCutMat_alpha "alpha" ///< Friction coefficient (DP model)
#define _IFT_DruckerPragerCutMat_alphapsi "alphapsi" /// Dilatancy coefficient (DP model)
#define _IFT_DruckerPragerCutMat_h "h" ///< Hardening modulus  (DP model)
#define _IFT_DruckerPragerCutMat_sigT "sigt" ///< Uniaxial tensile strength for cut-off, (Rankine plasticity model)
#define _IFT_DruckerPragerCutMat_omegaCrit "omega_crit" ///< Critical damage
#define _IFT_DruckerPragerCutMat_a "a" ///< Exponent in damage law
#define _IFT_DruckerPragerCutMat_yieldTol "yieldtol" ///< Tolerance of the error in the yield criterion
#define _IFT_DruckerPragerCutMat_newtonIter "newtoniter" ///< Maximum number of iterations in lambda search
#define _IFT_DruckerPragerCutMat_tau0 "tau0" ///< Initial yield stress under pure shear (DP model)
//@}

namespace oofem {
class GaussPoint;
class Domain;

/**
 * This class implements an isotropic elasto-plasto-damage material
 * with Drucker-Prager yield condition, tension cut-off, non-associated flow rule,
 * linear isotropic hardening and isotropic damage.
 */
class DruckerPragerCutMat : public MPlasticMaterial2
{
protected:
    // Reference to the basic elastic material.
    //LinearElasticMaterial *linearElasticMaterial;

    /// Elastic shear modulus.
    double G;

    /// Elastic bulk modulus.
    double K;

    /// Hardening modulus.
    double H;

    /// Uniaxial tensile strength for cut-off.
    double sigT;

    /// Initial yield stress under pure shear.
    double tau0;

    /// Friction coefficient.
    double alpha;

    ///Dilatancy coefficient (allowing non-associated plasticity).
    double alphaPsi;

    /// Tolerance of the error in the yield criterion.
    double yieldTol;

    /// Maximum number of iterations in lambda search.
    int newtonIter;

    /// Maximum damage value.
    double omegaCrit;

    /// Parameter for damage computation from cumulative plastic strain
    double a;

public:
    DruckerPragerCutMat(int n, Domain * d);
    virtual ~DruckerPragerCutMat();

    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual const char *giveClassName() const { return "DruckerPragerCutMat"; }
    virtual const char *giveInputRecordName() const { return _IFT_DruckerPragerCutMat_Name; }

    /// Returns a reference to the basic elastic material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual int giveSizeOfFullHardeningVarsVector() { return 4; }
    virtual int giveSizeOfReducedHardeningVarsVector(GaussPoint *) const { return 4; } //cummulative strain = one per each surface

protected:
    virtual int giveMaxNumberOfActiveYieldConds(GaussPoint *gp) { return 3; } //normally one less than number of all conditions

    virtual double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector, const FloatArray &strainSpaceHardeningVariables);

    virtual void computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector, const FloatArray &stressSpaceHardeningVars);

    /// Computes second derivative of yield/loading function with respect to stress
    virtual void computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int isurf, GaussPoint *gp, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables);

    virtual void computeReducedElasticModuli(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep);

    /// Functions related to hardening
    virtual int hasHardening() { return 1; }

    /// Compute dot(kappa_1), dot(kappa_2) etc.
    virtual void computeStrainHardeningVarsIncrement(FloatArray &answer, GaussPoint *gp, const FloatArray &stress, const FloatArray &dlambda, const FloatArray &dplasticStrain, const IntArray &activeConditionMap);

    /// Computes the derivative of yield/loading function with respect to kappa_1, kappa_2 etc.
    virtual void computeKGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables);

    /// computes mixed derivative of load function with respect to stress and hardening variables
    virtual void computeReducedSKGradientMatrix(FloatMatrix &gradientMatrix, int isurf, GaussPoint *gp, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables);

    /// computes dk(i)/dsig(j) gradient matrix
    virtual void computeReducedHardeningVarsSigmaGradient(FloatMatrix &answer, GaussPoint *gp, const IntArray &activeConditionMap, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVars, const FloatArray &dlambda);

    /// computes dKappa_i/dLambda_j
    virtual void computeReducedHardeningVarsLamGradient(FloatMatrix &answer, GaussPoint *gp, int actSurf, const IntArray &activeConditionMap, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVars, const FloatArray &dlambda);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
};
} // end namespace oofem
#endif // druckerpragercatmat_h
