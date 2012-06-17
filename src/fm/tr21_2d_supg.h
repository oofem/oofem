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

#ifndef tr21_2d_supg_h
#define tr21_2d_supg_h

#include "supgelement2.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "fei2dtrquad.h"
#include "fei2dtrlin.h"

#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "leplic.h"
#include "levelsetpcs.h"

namespace oofem {
/**
 * Class representing 2d triangular element  with quadratic velocity
 * and linear pressure approximation for solving incompressible fluid problems
 * with SUPG solver.
 */
class TR21_2D_SUPG : public SUPGElement2, public LevelSetPCSElementInterface, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI2dTrQuad velocityInterpolation;
    static FEI2dTrLin pressureInterpolation;
    IntArray pressureDofManArray;

public:
    TR21_2D_SUPG(int n, Domain *aDomain);
    virtual ~TR21_2D_SUPG();

    virtual FEInterpolation *giveInterpolation();
    virtual FEInterpolation *giveInterpolation(DofIDItem id);

    // definition
    virtual const char *giveClassName() const { return "TR21_2D_SUPG"; }
    virtual classType giveClassID() const { return TR21_2D_SUPGClass; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_2; }
    virtual MaterialMode giveMaterialMode() { return _2dFlow; }

    virtual void giveElementDofIDMask(EquationID, IntArray & answer) const;
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual int computeNumberOfDofs(EquationID ut);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void updateYourself(TimeStep *tStep);
    /// Used to check consistency and initialize some element geometry data (area,b,c).
    virtual int checkConsistency();

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual double LS_PCS_computeF(LevelSetPCS *ls, TimeStep *tStep);
    virtual void LS_PCS_computedN(FloatMatrix &answer);
    virtual double LS_PCS_computeVolume();
    virtual void LS_PCS_computeVolume(double &answer,  const FloatArray **coordinates);
    virtual double LS_PCS_computeS(LevelSetPCS *ls, TimeStep *tStep);
    virtual void LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi);


    virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }


    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }

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

    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *atTime);
    // Graphics output
    //void drawYourself(oofegGraphicContext&);
    virtual void drawRawGeometry(oofegGraphicContext &);
    virtual void drawScalar(oofegGraphicContext &context);
    //virtual void drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual double computeCriticalTimeStep(TimeStep *tStep);

    // three terms for computing their norms due to computing t_supg
    virtual void computeAdvectionTerm(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeAdvectionDeltaTerm(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeMassDeltaTerm(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeLSICTerm(FloatMatrix &answer, TimeStep *atTime);

    virtual Interface *giveInterface(InterfaceType);

protected:
    virtual void giveLocalVelocityDofMap (IntArray &map);
    virtual void giveLocalPressureDofMap (IntArray &map);

    virtual void computeGaussPoints();
    virtual void computeNuMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeUDotGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime);
    virtual void computeBMatrix(FloatMatrix &anwer, GaussPoint *gp);
    virtual void computeDivUMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeNpMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeGradPMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeDivTauMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime);
    virtual void computeGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime);
    virtual int  giveNumberOfSpatialDimensions();
    virtual double computeVolumeAround(GaussPoint *aGaussPoint);
    virtual void initGeometry();

    virtual void updateStabilizationCoeffs(TimeStep *tStep);

    virtual int giveTermIntergationRuleIndex(CharType termType);
};
} // end namespace oofem
#endif // tr21_2d_supg_h
