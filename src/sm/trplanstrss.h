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

#ifndef trplanstrss_h
#define trplanstrss_h

#include "nlstructuralelement.h"
#include "fei2dtrlin.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"
#include "directerrorindicatorrc.h"
#include "eleminterpmapperinterface.h"
#include "zzerrorestimator.h"
#include "mmashapefunctprojection.h"
#include "huertaerrorestimator.h"
#include "gausspoint.h"

#define _IFT_TrPlaneStress2d_Name "trplanestress2d"

namespace oofem {
/**
 * This class implements an triangular three-node plane-stress
 * elasticity finite element. Each node has 2 degrees of freedom.
 * Element has 3 nodes and 6 DoFs.
 * Tasks:
 * - calculating its B,D,N matrices and dV.
 */
class TrPlaneStress2d : public NLStructuralElement, public ZZNodalRecoveryModelInterface,
public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
public SpatialLocalizerInterface,
public EIPrimaryUnknownMapperInterface,
public ZZErrorEstimatorInterface, public MMAShapeFunctProjectionInterface,
public HuertaErrorEstimatorInterface
{
protected:
    static FEI2dTrLin interp;
    double area;

public:
    TrPlaneStress2d(int n, Domain * d);
    virtual ~TrPlaneStress2d() { }

    virtual FEInterpolation *giveInterpolation() const { return & interp; }
    virtual int computeNumberOfDofs() { return 6; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);
    virtual double giveCharacteristicSize(GaussPoint *gp, FloatArray &normalToCrackPlane, ElementCharSizeMethod method);
    virtual double giveParentElSize() const { return 0.5; }

    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual Interface *giveInterface(InterfaceType);

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawSpecial(oofegGraphicContext &gc, TimeStep *tStep);
#endif

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_TrPlaneStress2d_Name; }
    virtual const char *giveClassName() const { return "TrPlaneStress2d"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _PlaneStress; }

    //void ZZNodalRecoveryMI_computeNValProduct (FloatArray& answer, InternalStateType type, TimeStep* tStep);
    //void ZZNodalRecoveryMI_computeNNMatrix (FloatArray& answer, InternalStateType type);

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep);

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

    virtual void EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
                                                                       TimeStep *tStep, const FloatArray &lcoords,
                                                                       FloatArray &answer);

    // HuertaErrorEstimatorInterface
    virtual void HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                                  IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                                  HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                                  int &localNodeId, int &localElemId, int &localBcId,
                                                                  IntArray &controlNode, IntArray &controlDof,
                                                                  HuertaErrorEstimator :: AnalysisMode aMode);
    virtual void HuertaErrorEstimatorI_computeLocalCoords(FloatArray &answer, const FloatArray &coords)
    { computeLocalCoordinates(answer, coords); }
    virtual void HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
    { computeNmatrixAt(* ( gp->giveSubPatchCoordinates() ), answer); }

    virtual void MMAShapeFunctProjectionInterface_interpolateIntVarAt(FloatArray &answer, FloatArray &coords,
                                                                      coordType ct, nodalValContainerType &list,
                                                                      InternalStateType type, TimeStep *tStep);

protected:
    // edge load support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *gp);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);

    virtual void computeGaussPoints();

    virtual double giveArea();
    virtual FloatArray *GivebCoeff();
    virtual FloatArray *GivecCoeff();
    virtual int giveApproxOrder() { return 1; }
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 4; }
};
} // end namespace oofem
#endif // trplanstrss_h
