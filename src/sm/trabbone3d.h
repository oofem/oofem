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


#ifndef trabbone3d_h

#include "structuralmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "matconst.h"
#include "matstatus.h"
#include "strainvector.h"

#include "linearelasticmaterial.h"
#include "dictionr.h"

#include "structuralms.h"
#include "cltypes.h"

namespace oofem {

/**
 * This class implements associated Material Status to TrabBone3D.
 */
class TrabBone3DStatus : public StructuralMaterialStatus
{
protected:
    double kappa, tempKappa, dam, tempDam, tempPSED, tempTSED, tsed, beta;
    FloatArray tempPlasDef, plasDef, effectiveStress, tempEffectiveStress, plasFlowDirec;
    FloatMatrix smtrx, tangentMatrix, SSaTensor;

    /// Number of substeps in the last iteration.
    int nss;

public:
    TrabBone3DStatus(int n, Domain *d, GaussPoint *g);
    ~TrabBone3DStatus();

    void printOutputAt(FILE *file, TimeStep *tStep);

    double giveKappa();
    double giveTempKappa();
    double giveDam();
    double giveTempDam();
    double giveTempPSED();
    double giveTSED();
    double giveTempTSED();
    double giveBeta();
    int giveNsubsteps() {return nss;}

    const FloatArray *givePlasDef();
    const FloatArray *giveTempPlasDef();
    const FloatArray *giveTempEffectiveStress();
    const FloatArray *givePlasFlowDirec();
    const FloatMatrix *giveTangentMatrix();
    const FloatMatrix *giveSmtrx();
    const FloatMatrix *giveSSaTensor();

    void setTempKappa(double al) { tempKappa = al; }
    void setTempDam(double da) { tempDam = da; }
    void setTempPSED(double pse) { tempPSED = pse; }
    void setTempTSED(double tse) { tempTSED = tse; }
    void setBeta(double be) { beta = be; }
    void setTempEffectiveStress(FloatArray sc) { tempEffectiveStress = sc; }
    void setTempPlasDef(FloatArray epsip) { tempPlasDef = epsip; }
    void setPlasFlowDirec(FloatArray pfd) { plasFlowDirec = pfd; }
    void setSmtrx(FloatMatrix smt) { smtrx = smt; }
    void setTangentMatrix(FloatMatrix tmm) { tangentMatrix = tmm; }
    void setSSaTensor(FloatMatrix ssa) { SSaTensor = ssa; }
    void setNsubsteps(int n) { nss = n; }

    // definition
    const char *giveClassName() const { return "TrabBone3DStatus"; }
    classType giveClassID() const { return TrabBone3DStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};

/**
 * Trabecular bone material model in 3D.
 */
class TrabBone3D : public StructuralMaterial
{
protected:
    double m1, m2, rho, eps0, nu0, mu0, expk, expl, sig0Pos, sig0Neg, chi0Pos, chi0Neg, tau0, expq, expp;
    double plasHardFactor, expPlasHard, expDam, critDam, gamDens, tDens, JCrit;
    int printflag, abaqus, max_num_iter, max_num_substeps;
    double rel_yield_tol, strain_tol;

public:
    TrabBone3D(int n, Domain *d);

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }
    double evaluateCurrentYieldStress(const double kappa);
    double evaluateCurrentPlasticModulus(const double kappa);
    bool projectOnYieldSurface(double &tempKappa, FloatArray &tempEffectiveStress, FloatArray &tempPlasDef, const FloatArray &trialEffectiveStress, const FloatMatrix &elasticity, const FloatMatrix &compliance, TrabBone3DStatus *status);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);

    double computeDamageParam(double kappa, GaussPoint *gp);

    double computeDamage(GaussPoint *gp, TimeStep *atTime);

    /**
     * Computes the nonlocal cumulated plastic strain from its local form.
     * @param[out] kappa Return param, containing the corresponding cumulated plastic strain.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);

    void computePlasStrainEnerDensity(GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &totalStress);

    void computeDensificationStress(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *atTime);

    /// Constructs anisotropic compliance tensor.
    void constructAnisoComplTensor(FloatMatrix &answer);

    /// Constructs anisotropic fabric tensor.
    void constructAnisoFabricTensor(FloatMatrix &answer, const int posSignFlag);

    /// Construct tensor to adjust Norm.
    void constructNormAdjustTensor(FloatMatrix &answer);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode, GaussPoint *gp,
                                               TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual int hasMaterialModeCapability(MaterialMode);

    const char *giveClassName() const { return "TrabBone3D"; }
    classType giveClassID() const { return TrabBone3DClass; }

    IRResultType initializeFrom(InputRecord *ir);

    MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
};
} //end namespace oofem
#define trabbone3d_h
#endif
