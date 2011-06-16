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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   *******************************
//   *** CLASS Rankine material
//   *******************************

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
class GaussPoint;
class Domain;
class RankineMatStatus;

class RankineMat : public StructuralMaterial
{
    /*
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

protected:
    /// reference to the basic elastic material
    LinearElasticMaterial *linearElasticMaterial;

    /// Young's modulus
    double E;

    /// Poisson's ratio
    double nu;

    /// hardening modulus
    double H;

    /// initial (uniaxial) yield stress
    double sig0;

    /// parameter that controls damage evolution (a=0 turns damage off)
    double a;

public:
    RankineMat(int n, Domain *d);
    ~RankineMat() {; }
    double evalYieldFunction(const FloatArray &sigPrinc, const double kappa);
    double evalYieldStress(const double kappa);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, MaterialMode mode);
    double computeDamage(GaussPoint *gp, TimeStep *atTime);
    double computeDamageParam(double tempKappa);
    double computeDamageParamPrime(double tempKappa);
    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *atTime);
    /// specifies whether a given material mode is supported by this model
    int hasMaterialModeCapability(MaterialMode mode);

    /// reads the model parameters from the input file
    IRResultType initializeFrom(InputRecord *ir);

    // identification and auxiliary functions
    int hasNonLinearBehaviour()   { return 1; }
    const char *giveClassName() const { return "RankineMat"; }
    classType giveClassID()         const { return RankineMatClass; }

    /// returns a reference to the basic elastic material
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    /// determines whether the stiffness matrix is symmetric
    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return ( a == 0. ); }

    /// creates a new material status  corresponding to this class
    MaterialStatus *CreateStatus(GaussPoint *gp) const;

    /// evaluates the material stiffness matrix
    /*
     * void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
     *                                 MatResponseForm, MatResponseMode,
     *                                 GaussPoint * gp,
     *                                 TimeStep * atTime);
     */
    /// evaluates the stress
    void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);

protected:
    /// evaluates the consistent algorithmic stiffness under plane stress
    void givePlaneStressStiffMtrx(FloatMatrix & answer,
                                  MatResponseForm, MatResponseMode,
                                  GaussPoint * gp,
                                  TimeStep * atTime);
    /// executive method used by local and gradient version
    /// (with different parameters gprime)
    void evaluatePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                      MatResponseMode mode,
                                      GaussPoint *gp,
                                      TimeStep *atTime, double gprime);

    /// computes derivatives of final kappa with respect to final strain
    void computeEta(FloatArray &answer, RankineMatStatus *status);

    /**
     * Returns the integration point corresponding value in Reduced form.
     * @param answer contain corresponding ip value, zero sized if not available
     * @param aGaussPoint integration point
     * @param type determines the type of internal variable
     * @param type determines the type of internal variable
     * @returns nonzero if ok, zero if var not supported
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);

    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
     * @param type determines the internal variable requested (physical meaning)
     * @returns nonzero if ok or error is generated for unknown mat mode.
     */
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);

    /**
     * Returns the type of internal variable (scalar, vector, tensor,...).
     * @param type determines the type of internal variable
     * @returns type of internal variable
     */
    virtual InternalStateValueType giveIPValueType(InternalStateType type);

    /**
     * Returns the corresponding integration point  value size in Reduced form.
     * @param type determines the type of internal variable
     * @returns var size, zero if var not supported
     */
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
};

//=============================================================================


class RankineMatStatus : public StructuralMaterialStatus
{
protected:
    /// plastic strain (initial)
    FloatArray plasticStrain;

    /// plastic strain (final)
    FloatArray tempPlasticStrain;

    /// effective stress (initial)
    FloatArray effStress;

    /// effective stress (final)
    FloatArray tempEffStress;

    /// cumulative plastic strain (initial)
    double kappa;

    /// cumulative plastic strain (final)
    double tempKappa;

    /// increments of cumulative plastic strain
    /// associated with the first and secomnd principal stress
    /// (used in the case of vertex return, needed for stiffness)
    double dKappa1, dKappa2;

    /// damage (initial)
    double damage;

    /// damage (final)
    double tempDamage;

    /// tangent shear stiffness (needed for tangent matrix)
    double tanG;

#ifdef keep_track_of_dissipated_energy
    /// density of total work done by stresses on strain increments
    double stressWork;
    /// non-equilibrated density of total work done by stresses on strain increments
    double tempStressWork;
    /// density of dissipated work
    double dissWork;
    /// non-equilibrated density of dissipated work
    double tempDissWork;
#endif

public:
    RankineMatStatus(int n, Domain *d, GaussPoint *g);
    ~RankineMatStatus();

    void givePlasticStrain(FloatArray &answer)
    { answer = plasticStrain; }
    /*
     * void giveTrialStressDev(FloatArray &answer)
     * { answer = trialStressD; }
     *
     * void giveTrialStressVol(double &answer)
     * { answer = trialStressV; }
     */

    double giveDamage()
    { return damage; }

    double giveTempDamage()
    { return tempDamage; }

    double giveCumulativePlasticStrain()
    { return kappa; }

    double giveTempCumulativePlasticStrain()
    { return tempKappa; }

    double giveDKappa(int i)
    { if ( i == 1 ) { return dKappa1; } else { return dKappa2; } }

    double giveTangentShearStiffness()
    { return tanG; }

    void giveEffectiveStress(FloatArray &answer)
    { answer = effStress; }

    void giveTempEffectiveStress(FloatArray &answer)
    { answer = tempEffStress; }

    void letTempPlasticStrainBe(FloatArray values)
    { tempPlasticStrain = values; }

    void letEffectiveStressBe(FloatArray values)
    { effStress = values; }

    void letTempEffectiveStressBe(FloatArray values)
    { tempEffStress = values; }

    void setTempCumulativePlasticStrain(double value)
    { tempKappa = value; }

    void setDKappa(double val1, double val2)
    { dKappa1 = val1;
      dKappa2 = val2; }

    void setTempDamage(double value)
    { tempDamage = value; }

    void setTangentShearStiffness(double value)
    { tanG = value; }

    const FloatArray *givePlasDef()
    { return & plasticStrain; }

    /// prints the output variables into the *.out file
    void printOutputAt(FILE *file, TimeStep *tStep);

    /// initializes the temporary status
    virtual void initTempStatus();

    /// updates the state after a new equilibrium state has been reached
    virtual void updateYourself(TimeStep *);

    /// saves the current context(state) into a stream
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /// restores the state from a stream
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

#ifdef keep_track_of_dissipated_energy
    /// Returns the density of total work of stress on strain increments
    double giveStressWork() { return stressWork; }
    /// Returns the temp density of total work of stress on strain increments
    double giveTempStressWork() { return tempStressWork; }
    /// Sets the density of total work of stress on strain increments to given value
    void setTempStressWork(double w) { tempStressWork = w; }
    /// Returns the density of dissipated work
    double giveDissWork() { return dissWork; }
    /// Returns the density of temp dissipated work
    double giveTempDissWork() { return tempDissWork; }
    /// Sets the density of dissipated work to given value
    void setTempDissWork(double w) { tempDissWork = w; }
    /// computes the increment of total stress work and of dissipated work
    /// (gf is the dissipation density per unit volume at complete failure,
    /// it is needed only to determine which extremely small dissipation
    /// can be set to zero to get clean results, but parameter gf can be
    /// set to zero if not available)
    void computeWork(GaussPoint *, MaterialMode mode, double gf);
#endif

    /// identifies this class by its name
    const char *giveClassName() const { return "RankineMatStatus"; }

    /// identifies this class by its ID number
    classType             giveClassID() const
    { return RankineMatStatusClass; }
};
} // end namespace oofem
#endif // rankinemat_h
