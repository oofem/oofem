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

#ifndef tr21_2d_supg_h
#define tr21_2d_supg_h

#include "supgelement2.h"
#include "intarray.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "leplic.h"
#include "levelsetpcs.h"

#define _IFT_TR21_2D_SUPG_Name "tr21supg"

namespace oofem {
class FEI2dTrQuad;
class FEI2dTrLin;

/**
 * Class representing 2d triangular element with quadratic velocity
 * and linear pressure approximation for solving incompressible fluid problems
 * with SUPG solver.
 */
class TR21_2D_SUPG : public SUPGElement2, public LevelSetPCSElementInterface, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI2dTrQuad velocityInterpolation;
    static FEI2dTrLin pressureInterpolation;

public:
    TR21_2D_SUPG(int n, Domain * aDomain);

    FEInterpolation *giveInterpolation() const override;
    FEInterpolation *giveInterpolation(DofIDItem id) const override;

    // definition
    const char *giveClassName() const override { return "TR21_2D_SUPG"; }
    const char *giveInputRecordName() const override { return _IFT_TR21_2D_SUPG_Name; }
    MaterialMode giveMaterialMode() override { return _2dFlow; }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    int computeNumberOfDofs() override;
    void updateYourself(TimeStep *tStep) override;
    /// Used to check consistency and initialize some element geometry data (area,b,c).
    int checkConsistency() override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    double LS_PCS_computeF(LevelSetPCS *ls, TimeStep *tStep) override;
    void LS_PCS_computedN(FloatMatrix &answer) override;
    double LS_PCS_computeVolume() override;
    void LS_PCS_computeVolume(double &answer,  const FloatArray **coordinates);
    double LS_PCS_computeS(LevelSetPCS *ls, TimeStep *tStep) override;
    void LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi) override;


    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep) override;

    /// @name Helping functions for computing VOFFractions.
    //@{
    void computeIntersection(int iedge, FloatArray &intcoords, FloatArray &fi);
    void computeMiddlePointOnParabolicArc(FloatArray &answer, int iedge, FloatArray borderpoints);
    void computeCenterOf(FloatArray &C, FloatArray c, int dim);
    void computeQuadraticRoots(FloatArray Coeff, double &r1, double &r2);
    void computeCoordsOfEdge(FloatArray &answer, int iedge);
    void computeQuadraticFunct(FloatArray &answer, int iedge);
    void computeQuadraticFunct(FloatArray &answer, FloatArray line);
    //@{

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *tStep) override;
    // Graphics output
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
    //void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override {}
#endif

    double computeCriticalTimeStep(TimeStep *tStep) override;

    // three terms for computing their norms due to computing t_supg
    void computeAdvectionTerm(FloatMatrix &answer, TimeStep *tStep);
    void computeAdvectionDeltaTerm(FloatMatrix &answer, TimeStep *tStep);
    void computeMassDeltaTerm(FloatMatrix &answer, TimeStep *tStep);
    void computeLSICTerm(FloatMatrix &answer, TimeStep *tStep);

    Interface *giveInterface(InterfaceType) override;

protected:
    void giveLocalVelocityDofMap(IntArray &map) override;
    void giveLocalPressureDofMap(IntArray &map) override;

    void computeGaussPoints() override;
    void computeDeviatoricStress(FloatArray &answer, const FloatArray &eps, GaussPoint *gp, TimeStep *tStep) override;
    void computeTangent(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void computeNuMatrix(FloatMatrix &answer, GaussPoint *gp) override;
    void computeUDotGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep) override;
    void computeBMatrix(FloatMatrix &anwer, GaussPoint *gp) override;
    void computeDivUMatrix(FloatMatrix &answer, GaussPoint *gp) override;
    void computeNpMatrix(FloatMatrix &answer, GaussPoint *gp) override;
    void computeGradPMatrix(FloatMatrix &answer, GaussPoint *gp) override;
    void computeDivTauMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep) override;
    void computeGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep) override;
    int  giveNumberOfSpatialDimensions() override;
    double computeVolumeAround(GaussPoint *gp) override;
    void initGeometry();

    void updateStabilizationCoeffs(TimeStep *tStep) override;
};
} // end namespace oofem
#endif // tr21_2d_supg_h
