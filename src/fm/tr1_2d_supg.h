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

#ifndef tr1_2d_supg_h
#define tr1_2d_supg_h

#include "supgelement.h"
#include "primaryfield.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "leplic.h"
#include "levelsetpcs.h"

///@name Input fields for TR12DSUPG
//@{
#define _IFT_TR1_2D_SUPG_Name "tr1supg"
#define _IFT_Tr1SUPG_pvof "pvof"
#define _IFT_Tr1SUPG_vof "vof"
#define _IFT_Tr1SUPG2_mat0 "mat0"
#define _IFT_Tr1SUPG2_mat1 "mat1"
//@}

namespace oofem {
class FEI2dTrLin;

/**
 * Class representing 2d linear triangular element
 * for solving incompressible fluid with SUPG solver
 *
 * This class is similar to TR1_2D_SUPG2, but difference is in handling
 * multiple fluids. This class uses rule of mixture which interpolates the
 * properties using VOF value, requiring the use of twofluidmaterial class
 * as material model for this situation.
 */
class TR1_2D_SUPG : public SUPGElement,
public SpatialLocalizerInterface, public EIPrimaryFieldInterface,
public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
public LEPlicElementInterface, public LevelSetPCSElementInterface
{
protected:
    static FEI2dTrLin interp;

    //double a[3];
    double b [ 3 ];
    double c [ 3 ];
    double area = 0.;

public:
    TR1_2D_SUPG(int n, Domain * d);

    FEInterpolation *giveInterpolation() const override;

    void computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;
    void computeAdvectionTerm_MB(FloatArray &answer, TimeStep *tStep) override;
    void computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;
    void computeDiffusionTerm_MB(FloatArray &answer, TimeStep *tStep) override;
    void computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep) override;
    void computePressureTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;
    void computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;
    void computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;
    void computeAdvectionTerm_MC(FloatArray &answer, TimeStep *tStep) override;
    void computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;
    void computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep) override {
        answer.resize(3, 6);
        answer.zero();
    }
    void computeDiffusionTerm_MC(FloatArray &answer, TimeStep *tStep) override {
        answer.resize(3);
        answer.zero();
    }
    void computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;
    void computePressureTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;
    void computeBCRhsTerm_MB(FloatArray &answer, TimeStep *tStep) override;
    void computeBCRhsTerm_MC(FloatArray &answer, TimeStep *tStep) override;
    void computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep) override;

    void computeSlipWithFrictionBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *tStep) override;
    void computePenetrationWithResistanceBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *tStep) override;
    void computeOutFlowBCTerm_MB(FloatMatrix &answer, int side, TimeStep *tStep) override;

    void computeHomogenizedReinforceTerm_MB(FloatMatrix &answer, Load *load, TimeStep *tStep) override;
    void computeHomogenizedReinforceTerm_MC(FloatMatrix &answer, Load *load, TimeStep *tStep) override;

    void updateStabilizationCoeffs(TimeStep *tStep) override;
    double computeCriticalTimeStep(TimeStep *tStep) override;

    // definition
    const char *giveClassName() const override { return "TR1_2D_SUPG"; }
    const char *giveInputRecordName() const override { return _IFT_TR1_2D_SUPG_Name; }
    MaterialMode giveMaterialMode() override { return _2dFlow; }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    int computeNumberOfDofs() override;
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    void updateYourself(TimeStep *tStep) override;
    /// Used to check consistency and initialize some element geometry data (area,b,c).
    int checkConsistency() override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    Interface *giveInterface(InterfaceType) override;

    int  EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                               const FloatArray &coords, IntArray &dofId, ValueModeType mode,
                                               TimeStep *tStep) override;

    double computeLEPLICVolumeFraction(const FloatArray &n, const double p, LEPlic *matInterface, bool updFlag) override;
    void formMaterialVolumePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                const FloatArray &normal, const double p, bool updFlag) override;
    void formVolumeInterfacePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                 const FloatArray &normal, const double p, bool updFlag) override;
    double truncateMatVolume(const Polygon &matvolpoly, double &volume) override;
    void giveElementCenter(LEPlic *mat_interface, FloatArray &center, bool updFlag) override;
    void formMyVolumePoly(Polygon &myPoly, LEPlic *mat_interface, bool updFlag) override;
    Element *giveElement() override { return this; }
    double computeMyVolume(LEPlic *matInterface, bool updFlag) override;
    double computeVolumeAround(GaussPoint *gp) override;
    double computeCriticalLEPlicTimeStep(TimeStep *tStep) override;

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep) override;

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override;
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override;

    double LS_PCS_computeF(LevelSetPCS *ls, TimeStep *tStep) override;
    void LS_PCS_computedN(FloatMatrix &answer) override;
    double LS_PCS_computeVolume() override { return area; }
    double LS_PCS_computeS(LevelSetPCS *ls, TimeStep *tStep) override;
    void LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *tStep) override;
    // Graphics output
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
    //void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override {}
#endif

    void printOutputAt(FILE *file, TimeStep *tStep) override;

protected:
    void giveLocalVelocityDofMap(IntArray &map) override;
    void giveLocalPressureDofMap(IntArray &map) override;

    void computeNMtrx(FloatArray &answer, GaussPoint *gp);
    void computeGaussPoints() override;

    void computeDeviatoricStrain(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;
    void computeDeviatoricStress(FloatArray &answer, const FloatArray &eps, GaussPoint *gp, TimeStep *tStep) override;
    void computeTangent(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    virtual void initGeometry();
};
} // end namespace oofem
#endif // tr1_2d_supg_h
