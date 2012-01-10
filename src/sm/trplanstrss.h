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

#ifndef trplanstrss_h
#define trplanstrss_h

#include "nlstructuralelement.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"

#include "directerrorindicatorrc.h"
#include "eleminterpmapperinterface.h"
#include "zzerrorestimator.h"
#include "mmashapefunctprojection.h"
#include "huertaerrorestimator.h"

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
    public DirectErrorIndicatorRCInterface, public EIPrimaryUnknownMapperInterface,
    public ZZErrorEstimatorInterface, public ZZRemeshingCriteriaInterface, public MMAShapeFunctProjectionInterface,
    public HuertaErrorEstimatorInterface, public HuertaRemeshingCriteriaInterface
{
protected:
    double area;
    int numberOfGaussPoints;

public:
    TrPlaneStress2d(int n, Domain *d);
    ~TrPlaneStress2d() { }

    virtual int computeNumberOfDofs(EquationID ut) { return 6; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);
    double giveCharacteristicSize(GaussPoint *gp, FloatArray &normalToCrackPlane, ElementCharSizeMethod method);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }
    double computeVolumeAround(GaussPoint *gp);

    Interface *giveInterface(InterfaceType);

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void drawScalar(oofegGraphicContext &context);
    virtual void drawSpecial(oofegGraphicContext &);
    //void drawInternalState(oofegGraphicContext&);
#endif

    // definition & identification
    const char *giveClassName() const { return "TrPlaneStress2d"; }
    classType giveClassID() const { return TrPlaneStress2dClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }

    //void ZZNodalRecoveryMI_computeNValProduct (FloatArray& answer, InternalStateType type, TimeStep* tStep);
    //void ZZNodalRecoveryMI_computeNNMatrix (FloatArray& answer, InternalStateType type);
    int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    Element *ZZNodalRecoveryMI_giveElement() { return this; }
    void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint,
                                                             InternalStateType type);

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    int SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
    int SPRNodalRecoveryMI_giveNumberOfIP();
    //void SPRNodalRecoveryMI_giveIPValue (FloatArray& answer, int ipNum, InternalStateType type);
    void SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp);
    SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

    virtual double DirectErrorIndicatorRCI_giveCharacteristicSize();
    virtual int EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType u,
                                                                 TimeStep *stepN, const FloatArray &coords,
                                                                 FloatArray &answer);
    virtual void EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer);

    // ZZErrorEstimatorInterface
    virtual Element *ZZErrorEstimatorI_giveElement() { return this; }
    virtual void ZZErrorEstimatorI_computeEstimatedStressInterpolationMtrx(FloatMatrix &answer, GaussPoint *gp,
                                                                           InternalStateType type)
    { ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(answer, gp, type); }

    // HuertaErrorEstimatorInterface
    virtual Element *HuertaErrorEstimatorI_giveElement() { return this; }

    virtual void HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                                  IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                                  HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                                  int &localNodeId, int &localElemId, int &localBcId,
                                                                  IntArray &controlNode, IntArray &controlDof,
                                                                  HuertaErrorEstimator :: AnalysisMode aMode);
    void HuertaErrorEstimatorI_computeLocalCoords(FloatArray &answer, const FloatArray &coords)
    { computeLocalCoordinates(answer, coords); }
    void HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
    { computeNmatrixAt(aGaussPoint, answer); }

    // ZZRemeshingCriteriaInterface
    virtual double ZZRemeshingCriteriaI_giveCharacteristicSize() { return DirectErrorIndicatorRCI_giveCharacteristicSize(); }
    virtual int ZZRemeshingCriteriaI_givePolynOrder() { return 1; };

    // HuertaRemeshingCriteriaInterface
    virtual double HuertaRemeshingCriteriaI_giveCharacteristicSize() { return DirectErrorIndicatorRCI_giveCharacteristicSize(); };
    virtual int HuertaRemeshingCriteriaI_givePolynOrder() { return 1; };

    virtual void MMAShapeFunctProjectionInterface_interpolateIntVarAt(FloatArray &answer, FloatArray &coords,
                                                                      coordType ct, nodalValContainerType &list,
                                                                      InternalStateType type, TimeStep *tStep);

    integrationDomain giveIntegrationDomain() { return _Triangle; }
    MaterialMode giveMaterialMode() { return _PlaneStress; }

protected:
    // edge load support
    void computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *gp);
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    void computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *gp, int);

    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeGaussPoints();

    virtual double giveArea();
    virtual FloatArray *GivebCoeff();
    virtual FloatArray *GivecCoeff();
    int giveApproxOrder() { return 1; }
    /**
     * Return desired number of integration points for consistent mass matrix
     * computation, if required.
     */
    int giveNumberOfIPForMassMtrxIntegration() { return 4; }
};
} // end namespace oofem
#endif // trplanstrss_h
