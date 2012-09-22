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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef druckerpragercatmat_h
#define druckerpragercatmat_h

#include "structuralmaterial.h"
#include "structuralms.h"
#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
class GaussPoint;
class Domain;

/**
 * This class implements an isotropic elasto-plasto-damage material
 * with Drucker-Prager yield condition, tension cut-off, non-associated flow rule,
 * linear isotropic hardening and isotropic damage.
 */
class DruckerPragerCutMat : public StructuralMaterial
{
protected:
    /// Reference to the basic elastic material.
    LinearElasticMaterial *linearElasticMaterial;

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
    DruckerPragerCutMat(int n, Domain *d);
    virtual ~DruckerPragerCutMat();

    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, MaterialMode mode);
    double computeDamage(GaussPoint *gp, TimeStep *atTime);
    double computeDamageParam(double tempKappa);
    double computeDamageParamPrime(double tempKappa);
    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *atTime);

    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual const char *giveClassName() const { return "DruckerPragerCutMat"; }
    virtual classType giveClassID() const { return DruckerPragerCutMatClass; }

    /// Returns a reference to the basic elastic material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                       MatResponseForm form, MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep);

protected:
    /// Evaluates the stress from Green-Lagrange strain E.
    void giveRealStressVectorComputedFromStrain(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                                const FloatArray &E, TimeStep *tStep);

    void give3dSSMaterialStiffnessMatrix(FloatMatrix &answer,
                                         MatResponseForm form, MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *tStep);
    virtual void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);

    virtual InternalStateValueType giveIPValueType(InternalStateType type);

    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
};

//=============================================================================


class DruckerPragerCutMatStatus : public StructuralMaterialStatus
{
protected:
    /// Plastic strain (initial).
    FloatArray plasticStrain;

    /// Plastic strain (temporary and final after iterations).
    FloatArray tempPlasticStrain;

    /// Deviatoric trial stress - needed for tangent stiffness.
    FloatArray trialStressD;

    double trialStressV;

    FloatArray effStress;
    FloatArray tempEffStress;

    /// Cumulative plastic strain (initial).
    double kappa;

    /// Cumulative plastic strain (final).
    double tempKappa;

    double tempDamage, damage;
    
    /// Active surface (None, DruckerPrager, Cut-off)
    int activeSurface;

public:
    DruckerPragerCutMatStatus(int n, Domain *d, GaussPoint *g);
    virtual ~DruckerPragerCutMatStatus();

    void givePlasticStrain(FloatArray &answer) { answer = plasticStrain; }

    void giveTrialStressDev(FloatArray &answer) { answer = trialStressD; }

    void giveTrialStressVol(double &answer) { answer = trialStressV; }
    double giveDamage() { return damage; }
    double giveTempDamage() { return tempDamage; }

    double giveCumulativePlasticStrain() { return kappa; }
    double giveTempCumulativePlasticStrain() { return tempKappa; }
    int giveActiveSurface() { return activeSurface; }

    void giveTempEffectiveStress(FloatArray &answer) { answer = tempEffStress; }
    void giveEffectiveStress(FloatArray &answer) { answer = effStress; }

    void letTempPlasticStrainBe(FloatArray &values) { tempPlasticStrain = values; }

    void letTrialStressDevBe(FloatArray &values) { trialStressD = values; }

    void letEffectiveStressBe(FloatArray &values) { effStress = values; }

    void letTempEffectiveStressBe(FloatArray &values) { tempEffStress = values; }

    void setTrialStressVol(double value) { trialStressV = value; }

    void setTempCumulativePlasticStrain(double value) { tempKappa = value; }

    void setTempDamage(double value) { tempDamage = value; }
    
    void setActiveSurface(int as) { activeSurface = as; }

    const FloatArray *givePlasDef() { return & plasticStrain; }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "DruckerPragerCutMatStatus"; }
    virtual classType giveClassID() const { return DruckerPragerCutMatStatusClass; }
};
} // end namespace oofem
#endif // druckerpragercatmat_h
