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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef rankinemat_h
#define rankinemat_h

#include "structuralmaterial.h"
#include "structuralms.h"
#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

// this turns on or off a bunch of internal variables
// that allow tracking the distribution of dissipated energy
// (can be turned off if such information is not needed)
#define keep_track_of_dissipated_energy
//#undefine keep_track_of_dissipated_energy

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
    LinearElasticMaterial *linearElasticMaterial;

    /// Young's modulus.
    double E;

    /// Poisson's ratio.
    double nu;

    /// Hardening modulus.
    double H;

    /// Initial (uniaxial) yield stress.
    double sig0;

    /// Parameter that controls damage evolution (a=0 turns damage off).
    double a;

public:
    RankineMat(int n, Domain *d);
    virtual ~RankineMat() { }

    double evalYieldFunction(const FloatArray &sigPrinc, const double kappa);
    double evalYieldStress(const double kappa);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, MaterialMode mode);
    double computeDamage(GaussPoint *gp, TimeStep *atTime);
    double computeDamageParam(double tempKappa);
    double computeDamageParamPrime(double tempKappa);
    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *atTime);

    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual IRResultType initializeFrom(InputRecord *ir);

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour() { return 1; }
    virtual const char *giveClassName() const { return "RankineMat"; }
    virtual classType giveClassID() const { return RankineMatClass; }

    /// Returns a reference to the basic elastic material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return ( a == 0. ); }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                              const FloatArray &reducesStrain, TimeStep *tStep);

protected:
    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer,
                                  MatResponseForm form, MatResponseMode mode,
                                  GaussPoint *gp,
                                  TimeStep *tStep);
    /**
     * Executive method used by local and gradient version.
     * (with different parameters gprime)
     */
    void evaluatePlaneStressStiffMtrx(FloatMatrix &answer,
                                      MatResponseForm form, MatResponseMode mode,
                                      GaussPoint *gp,
                                      TimeStep *tStep, double gprime);

    /// Computes derivatives of final kappa with respect to final strain.
    void computeEta(FloatArray &answer, RankineMatStatus *status);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *gp);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);

    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
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
    double kappa;

    /// Cumulative plastic strain (final).
    double tempKappa;

    /**
     * Increments of cumulative plastic strain
     * associated with the first and secomnd principal stress
     * (used in the case of vertex return, needed for stiffness)
     */
    double dKappa1, dKappa2;

    /// Damage (initial).
    double damage;

    /// Damage (final).
    double tempDamage;

    /// Tangent shear stiffness (needed for tangent matrix).
    double tanG;

#ifdef keep_track_of_dissipated_energy
    /// Density of total work done by stresses on strain increments.
    double stressWork;
    /// Non-equilibrated density of total work done by stresses on strain increments.
    double tempStressWork;
    /// Density of dissipated work.
    double dissWork;
    /// Non-equilibrated density of dissipated work.
    double tempDissWork;
#endif

public:
    RankineMatStatus(int n, Domain *d, GaussPoint *g);
    virtual ~RankineMatStatus();

    void givePlasticStrain(FloatArray &answer) { answer = plasticStrain; }

    double giveDamage() { return damage; }
    double giveTempDamage() { return tempDamage; }

    double giveCumulativePlasticStrain() { return kappa; }
    double giveTempCumulativePlasticStrain() { return tempKappa; }

    double giveDKappa(int i)
    { if ( i == 1 ) { return dKappa1; } else { return dKappa2; } }

    double giveTangentShearStiffness()
    { return tanG; }

    void giveEffectiveStress(FloatArray &answer) { answer = effStress; }
    void giveTempEffectiveStress(FloatArray &answer) { answer = tempEffStress; }

    void letTempPlasticStrainBe(FloatArray &values) { tempPlasticStrain = values; }

    void letEffectiveStressBe(FloatArray &values) { effStress = values; }

    void letTempEffectiveStressBe(FloatArray &values) { tempEffStress = values; }

    void setTempCumulativePlasticStrain(double value) { tempKappa = value; }

    void setDKappa(double val1, double val2) {
        dKappa1 = val1;
        dKappa2 = val2;
    }

    void setTempDamage(double value) { tempDamage = value; }

    void setTangentShearStiffness(double value) { tanG = value; }

    const FloatArray *givePlasDef() { return & plasticStrain; }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

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
    void computeWork(GaussPoint *gp, MaterialMode mode, double gf);
#endif

    virtual const char *giveClassName() const { return "RankineMatStatus"; }
    virtual classType giveClassID() const { return RankineMatStatusClass; }
};
} // end namespace oofem
#endif // rankinemat_h
