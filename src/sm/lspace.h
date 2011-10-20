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

#ifndef lspace_h
#define lspace_h

#include "nlstructuralelement.h"
#include "fei3dhexalin.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "spatiallocalizer.h"
#include "eleminterpmapperinterface.h"
#include "huertaerrorestimator.h"

namespace oofem {
/**
 * This class implements a Linear 3d 8-node finite element for stress analysis.
 * Each node has 3 degrees of freedom.
 *
 * One single additional attribute is needed for Gauss integration purpose :
 * 'jacobianMatrix'. This 3x3 matrix contains polynomials.
 * TASKS :
 * - Calculating its Gauss points.
 * - Calculating its B,D,N matrices and dV.
 */
class LSpace  : public NLStructuralElement, public ZZNodalRecoveryModelInterface,
    public SPRNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface,
    public SpatialLocalizerInterface,
    public EIPrimaryUnknownMapperInterface,
    public HuertaErrorEstimatorInterface, public HuertaRemeshingCriteriaInterface
{
protected:
    static FEI3dHexaLin interpolation;
    int numberOfGaussPoints;

public:
    LSpace(int n, Domain *d);
    ~LSpace() { }

    virtual int computeNumberOfDofs(EquationID ut) { return 24; }
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    double computeVolumeAround(GaussPoint *gp);
    FEInterpolation *giveInterpolation() { return & interpolation; }

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);

    Interface *giveInterface(InterfaceType it);
    virtual int testElementExtension(ElementExtension ext)
    { return ( ( ( ext == Element_EdgeLoadSupport ) || ( ext == Element_SurfaceLoadSupport ) ) ? 1 : 0 ); }

    Element *ZZNodalRecoveryMI_giveElement() { return this; }

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    int SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
    int SPRNodalRecoveryMI_giveNumberOfIP();
    //void SPRNodalRecoveryMI_giveIPValue (FloatArray& answer, int ipNum, InternalStateType type);
    void SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp);
    SPRPatchType SPRNodalRecoveryMI_givePatchType();

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }

    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

    virtual int EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType u,
                                                                 TimeStep *stepN, const FloatArray &coords,
                                                                 FloatArray &answer);
    virtual void EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer);

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

    // HuertaRemeshingCriteriaInterface
    virtual double HuertaRemeshingCriteriaI_giveCharacteristicSize();
    virtual int HuertaRemeshingCriteriaI_givePolynOrder() { return 1; };

    // definition & identification
    const char *giveClassName() const { return "LSpace"; }
    classType giveClassID() const { return LSpaceClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_hexa_1; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType type);
    virtual void drawScalar(oofegGraphicContext &context);
    virtual void drawSpecial(oofegGraphicContext &);
    void drawTriad(FloatArray &, int isurf);
#endif

    integrationDomain giveIntegrationDomain() { return _Cube; }
    MaterialMode giveMaterialMode();

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    void computeBFmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *gp, int i);
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeGaussPoints();
    int giveApproxOrder() { return 1; }
    int giveNumberOfIPForMassMtrxIntegration() { return 8; }

    /**
     * @name Edge load support
     */
    //@{
    void computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
    //@}

    /**
     * @name Surface load support
     */
    //@{
    virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    virtual void giveSurfaceDofMapping(IntArray &answer, int iSurf) const;
    virtual IntegrationRule *GetSurfaceIntegrationRule(int iSurf);
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf);
    virtual void computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf);
    virtual int computeLoadLSToLRotationMatrix(FloatMatrix &answer, int iSurf, GaussPoint *gp);
    //@}
};
} // end namespace oofem
#endif // lspace_h
