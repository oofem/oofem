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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "matconst.h"
#include "matstatus.h"
#include "sm/Materials/structuralms.h"
#include "cltypes.h"
#include "floatmatrixf.h"

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


#define _IFT_TrabBone3D_gMin "gmin"
#define _IFT_TrabBone3D_gammaL "gammal"
#define _IFT_TrabBone3D_gammaP "gammap"
#define _IFT_TrabBone3D_tDens "tdens"
#define _IFT_TrabBone3D_densCrit "denscrit"

#define _IFT_TrabBone3D_printflag "printflag"
#define _IFT_TrabBone3D_max_num_iter "max_num_iter"
#define _IFT_TrabBone3D_max_num_substeps "max_num_substeps"
#define _IFT_TrabBone3D_rel_yield_tol "rel_yield_tol"
#define _IFT_TrabBone3D_strain_tol "strain_tol"
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
    double kappa = 0., tempKappa = 0;
    double dam = 0., tempDam = 0.;
    double tempPSED = 0., tempTSED = 0.;
    double tsed = 0.;
    double beta = 0.;

    FloatArrayF<6> tempPlasDef, plasDef;
    FloatArrayF<6> effectiveStress, tempEffectiveStress;
    FloatArrayF<6> plasFlowDirec, tempStrain;

    FloatMatrixF<6,6> SSaTensor;

    /// Densificator criterion
    double densG = 1.;

public:
    TrabBone3DStatus(GaussPoint *g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    
    double giveKappa() const { return kappa; }
    double giveTempKappa() const { return tempKappa; }
    double giveDam() const { return dam; }
    double giveTempDam() const { return tempDam; }
    double giveTempPSED() const { return tempPSED; }
    double giveTSED() const { return tsed; }
    double giveTempTSED() const { return tempTSED; }
    double giveBeta() const { return beta; }
    double giveDensG() const { return densG; }

    const FloatArrayF<6> &givePlasDef() const { return plasDef; }
    const FloatArrayF<6> &giveTempPlasDef() const { return tempPlasDef; }
    const FloatArrayF<6> &givePlasFlowDirec() const { return plasFlowDirec; }
    const FloatArrayF<6> &giveTempEffectiveStress() const { return tempEffectiveStress; }
    const FloatMatrixF<6,6> &giveSSaTensor() const { return SSaTensor; }

    void setTempKappa(double al) { tempKappa = al; }
    void setKappa(double values) { kappa = values; }
    void setTempDam(double da) { tempDam = da; }
    void setTempPSED(double pse) { tempPSED = pse; }
    void setTempTSED(double tse) { tempTSED = tse; }
    void setBeta(double be) { beta = be; }
    void setTempEffectiveStress(const FloatArrayF<6> &sc) { tempEffectiveStress = sc; }
    void setTempPlasDef(const FloatArrayF<6> &epsip) { tempPlasDef = epsip; }
    void setPlasFlowDirec(const FloatArrayF<6> &pfd) { plasFlowDirec = pfd; }
    void setSSaTensor(const FloatMatrixF<6,6> &ssa) { SSaTensor = ssa; }

    void setDensG(double g) { densG = g; }

    const char *giveClassName() const override { return "TrabBone3DStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/////////////////////////////////////////////////////////////////
////////////////TRABECULAR BONE 3D///////////////////////////////
/////////////////////////////////////////////////////////////////


class TrabBone3D : public StructuralMaterial
{
protected:
    double m1 = 0., m2 = 0.;
    double rho = 0.;
    double eps0 = 0., nu0 = 0., mu0 = 0.;
    double expk = 0., expl = 0.;
    double sig0Pos = 0., sig0Neg = 0.;
    double chi0Pos = 0., chi0 = 0., chi0Neg = 0.;
    double tau0 = 0.;
    double expq = 0., expp = 0.;
    double plasHardFactor = 0., expPlasHard = 0., expDam = 0., critDam = 0.;
    int printflag = 0, max_num_iter = 0;
    double rel_yield_tol = 0., strain_tol = 0.;
    /// Local coordinate system
    double x1 = 1., x2 = 0., x3 = 0.;
    double y1 = 0., y2 = 1., y3 = 0.;
    double z1 = 0., z2 = 0., z3 = 1.;
    /// Densificator properties
    double gammaL0 = 0., gammaP0 = 0., tDens = 0., densCrit = 0., rL = 0., rP = 0., gammaL = 0., gammaP = 0.;
    /// Viscosity parameter
    double viscosity = 0.;

public:
    TrabBone3D(int n, Domain *d);

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }
    double evaluateCurrentYieldStress(double kappa) const;
    double evaluateCurrentPlasticModulus(double kappa) const;
    double evaluateCurrentViscousStress(double deltaKappa, TimeStep *tStep) const;
    double evaluateCurrentViscousModulus(double deltaKappa, TimeStep *tStep) const;

    bool projectOnYieldSurface(double &tempKappa, FloatArrayF<6> &tempEffectiveStress, FloatArrayF<6> &tempPlasDef,
                               const FloatArrayF<6> &trialEffectiveStress,
                               const FloatMatrixF<6,6> &elasticity, const FloatMatrixF<6,6> &compliance,
                               TrabBone3DStatus *status, TimeStep *tStep, GaussPoint *gp, int lineSearchFlag) const;

    void performPlasticityReturn(GaussPoint *gp, const FloatArrayF<6> &strain, TimeStep *tStep) const;

    std::pair<FloatArrayF<6>, double> constructPlasFlowDirec(FloatMatrixF<6,6> &fabric, const FloatArrayF<6> &F, const FloatArrayF<6> &S) const;
    FloatMatrixF<6,6> constructDerivativeOfPlasFlowDirec(const FloatMatrixF<6,6> &fabric, const FloatArrayF<6> &F, const FloatArrayF<6> &S) const;
    double evaluatePlasCriterion(const FloatMatrixF<6,6> &fabric, const FloatArrayF<6> &F, const FloatArrayF<6> &stress, double kappa, double deltaKappa, TimeStep *tStep) const;

    double computeDamageParam(double kappa) const;
    double computeDamageParamPrime(double kappa) const;

    double computeDamage(GaussPoint *gp, TimeStep *tStep) const;

    virtual double computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const;

    void computePlasStrainEnerDensity(GaussPoint *gp, const FloatArrayF<6> &strain, const FloatArrayF<6> &stress) const;

    FloatArrayF<6> computeDensificationStress(GaussPoint *gp, const FloatArrayF<6> &totalStrain, TimeStep *tStep) const;

    /// Construct anisotropic compliance tensor.
    FloatMatrixF<6,6> constructAnisoComplTensor() const;
    /// Construct anisotropic stiffness tensor.
    FloatMatrixF<6,6> constructAnisoStiffnessTensor() const;

    /// Construct anisotropic fabric tensor.
    FloatMatrixF<6,6> constructAnisoFabricTensor() const;
    FloatArrayF<6> constructAnisoFtensor() const;

    FloatMatrixF<6,6> constructStiffnessTransformationMatrix() const;
    FloatMatrixF<6,6> constructFabricTransformationMatrix() const;
    /// Construct Tensor to adjust Norm.
    FloatMatrixF<6,6> constructNormAdjustTensor() const;

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode, GaussPoint * gp, TimeStep * tStep) const override;

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;

    const char *giveInputRecordName() const override { return _IFT_TrabBone3D_Name; }
    const char *giveClassName() const override { return "TrabBone3D"; }

    void initializeFrom(InputRecord &ir) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    double predictRelativeComputationalCost(GaussPoint *gp) override;
    double predictRelativeRedistributionCost(GaussPoint *gp) override;
};
} //end namespace oofem
#endif
