/*
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

#ifndef trabbone3d_h
#define trabbone3d_h

#include "../sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "matconst.h"
#include "matstatus.h"
#include "../sm/Materials/structuralms.h"
#include "cltypes.h"

///@name Input fields for TrabBone3D
//@{
#define _IFT_TrabBone3D_Name "trabbone3d"
#define _IFT_TrabBone3D_eps0 "eps0"
#define _IFT_TrabBone3D_nu0 "nu0"
#define _IFT_TrabBone3D_mu0 "mu0"
#define _IFT_TrabBone3D_expk "expk"
#define _IFT_TrabBone3D_expl "expl"

#define _IFT_TrabBone3D_m1 "m1"
#define _IFT_TrabBone3D_m2 "m2"
#define _IFT_TrabBone3D_rho "rho"

#define _IFT_TrabBone3D_sig0Pos "sig0pos"
#define _IFT_TrabBone3D_sig0Neg "sig0neg"
#define _IFT_TrabBone3D_chi0Pos "chi0pos"
#define _IFT_TrabBone3D_chi0Neg "chi0neg"
#define _IFT_TrabBone3D_tau0 "tau0"
#define _IFT_TrabBone3D_expp "expp"
#define _IFT_TrabBone3D_expq "expq"
#define _IFT_TrabBone3D_plasHardFactor "plashardfactor"
#define _IFT_TrabBone3D_expPlasHard "expplashard"

#define _IFT_TrabBone3D_expDam "expdam"
#define _IFT_TrabBone3D_critDam "critdam"

#define _IFT_TrabBone3D_x1 "x1"
#define _IFT_TrabBone3D_x2 "x2"
#define _IFT_TrabBone3D_x3 "x3"
#define _IFT_TrabBone3D_y1 "y1"
#define _IFT_TrabBone3D_y2 "y2"
#define _IFT_TrabBone3D_y3 "y3"
#define _IFT_TrabBone3D_viscosity "viscosity"
#define _IFT_TrabBone3D_yR "yr"
#define _IFT_TrabBone3D_kappaMax "kappamax"
#define _IFT_TrabBone3D_kappaMin "kappamin"
#define _IFT_TrabBone3D_kappaSlope "kappaslope"
#define _IFT_TrabBone3D_N "n"
#define _IFT_TrabBone3D_gMin "gmin"
#define _IFT_TrabBone3D_formulation "formulation"
#define _IFT_TrabBone3D_gammaL "gammal"
#define _IFT_TrabBone3D_gammaP "gammap"
#define _IFT_TrabBone3D_tDens "tdens"
#define _IFT_TrabBone3D_densCrit "denscrit"
#define _IFT_TrabBone3D_printflag "printflag"
#define _IFT_TrabBone3D_max_num_iter "max_num_iter"
#define _IFT_TrabBone3D_max_num_substeps "max_num_substeps"
#define _IFT_TrabBone3D_rel_yield_tol "rel_yield_tol"
#define _IFT_TrabBone3D_strain_tol "strain_tol"
#define _IFT_TrabBone3D_abaqus "abaqus"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to TrabBone3D (trabecular bone material).
 * It is attribute of matStatusDictionary at every GaussPoint, for which this material
 * is active.
 */
class TrabBone3DStatus : public StructuralMaterialStatus
{
protected:
    double kappa, tempKappa, dam, tempDam, tempPSED, tempTSED, tsed, beta;
    FloatArray tempPlasDef, plasDef, effectiveStress, tempEffectiveStress, plasFlowDirec, tempStrain;
    ;
    FloatMatrix smtrx, tangentMatrix, SSaTensor;
    /// Number of substeps in the last iteration.
    int nss;
    /// Densificator criterion
    double densG;


public:
    TrabBone3DStatus(int n, Domain * d, GaussPoint * g);

    virtual ~TrabBone3DStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    double giveKappa();
    double giveTempKappa();
    double giveDam();
    double giveTempDam();
    double giveTempPSED();
    double giveTSED();
    double giveTempTSED();
    double giveBeta();
    int giveNsubsteps() { return nss; }
    double giveDensG() { return densG; }

    const FloatArray &givePlasDef() const;
    const FloatArray &giveTempPlasDef() const;
    const FloatArray &giveTempEffectiveStress() const;
    const FloatArray &givePlasFlowDirec() const;
    const FloatMatrix &giveTangentMatrix() const;
    const FloatMatrix &giveSmtrx() const;
    const FloatMatrix &giveSSaTensor() const;

    void setTempKappa(double al) { tempKappa = al; }
    void setKappa(double values) { kappa = values; }
    void setTempDam(double da) { tempDam = da; }
    void setTempPSED(double pse) { tempPSED = pse; }
    void setTempTSED(double tse) { tempTSED = tse; }
    void setBeta(double be) { beta = be; }
    void setTempEffectiveStress(FloatArray &sc) { tempEffectiveStress = sc; }
    void setTempPlasDef(FloatArray &epsip) { tempPlasDef = epsip; }
    void setPlasFlowDirec(FloatArray &pfd) { plasFlowDirec = pfd; }
    void setSmtrx(FloatMatrix &smt) { smtrx = smt; }
    void setTangentMatrix(FloatMatrix &tmm) { tangentMatrix = tmm; }
    void setSSaTensor(FloatMatrix &ssa) { SSaTensor = ssa; }
    void setNsubsteps(int n)  { nss = n; }

    void setDensG(double g) { densG = g; }

    virtual const char *giveClassName() const { return "TrabBone3DStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};


/////////////////////////////////////////////////////////////////
////////////////TRABECULAR BONE 3D///////////////////////////////
/////////////////////////////////////////////////////////////////


class TrabBone3D : public StructuralMaterial
{
protected:
    double m1, m2, rho, eps0, nu0, mu0, expk, expl, sig0Pos, sig0Neg, chi0Pos, chi0, chi0Neg, tau0, expq, expp;
    double plasHardFactor, expPlasHard, expDam, critDam, pR;
    int printflag, abaqus, max_num_iter, max_num_substeps;
    double rel_yield_tol, strain_tol;
    /// Local coordinate system
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    /// Densificator properties
    double gammaL0, gammaP0, tDens, densCrit, rL, rP, gammaL, gammaP;
    /// Viscosity parameter
    double viscosity;
    /// Hadi post-yield function
    double yR, kappaMax, kappaMin, kappaSlope, N, gMin, formulation;
    double hardFactor;

public:
    TrabBone3D(int n, Domain * d);

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }
    double evaluateCurrentYieldStress(const double kappa);
    double evaluateCurrentPlasticModulus(const double kappa);
    double evaluateCurrentViscousStress(const double deltaKappa, TimeStep *tStep);
    double evaluateCurrentViscousModulus(const double deltaKappa, TimeStep *tStep);

    bool projectOnYieldSurface(double &tempKappa, FloatArray &tempEffectiveStress, FloatArray &tempPlasDef, const FloatArray &trialEffectiveStress, const FloatMatrix &elasticity, const FloatMatrix &compliance, TrabBone3DStatus *status, TimeStep *tStep, GaussPoint *gp, int lineSearchFlag);

    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep);

    void  constructPlasFlowDirec(FloatArray &answer, double &norm, FloatMatrix &fabric, FloatArray &F, FloatArray &S);
    void  constructDerivativeOfPlasFlowDirec(FloatMatrix &answer, FloatMatrix &fabric, FloatArray &F, FloatArray &S);
    double evaluatePlasCriterion(FloatMatrix &fabric, FloatArray &F, FloatArray &stress, double kappa, double deltaKappa, TimeStep *tStep);

    double computeDamageParam(double kappa);
    double computeDamageParamPrime(double kappa);

    double computeDamage(GaussPoint *gp, TimeStep *tStep);

    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);

    void computePlasStrainEnerDensity(GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &totalStress);

    void computeDensificationStress(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep);

    /// Construct anisotropic compliance tensor.
    void constructAnisoComplTensor(FloatMatrix &answer);
    /// Construct anisotropic stiffness tensor.
    void constructAnisoStiffnessTensor(FloatMatrix &answer);

    /// Construct anisotropic fabric tensor.
    void constructAnisoFabricTensor(FloatMatrix &answer);
    void constructAnisoFtensor(FloatArray &answer);

    void constructStiffnessTransformationMatrix(FloatMatrix &answer);
    void constructFabricTransformationMatrix(FloatMatrix &answer);
    /// Construct Tensor to adjust Norm.
    void constructNormAdjustTensor(FloatMatrix &answer);


    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode, GaussPoint *gp,
                                               TimeStep *tStep);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *,
                                         const FloatArray &, TimeStep *);

    virtual const char *giveInputRecordName() const { return _IFT_TrabBone3D_Name; }
    virtual const char *giveClassName() const { return "TrabBone3D"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual double predictRelativeComputationalCost(GaussPoint *gp);
    virtual double predictRelativeRedistributionCost(GaussPoint *gp);
};
} //end namespace oofem
#endif
