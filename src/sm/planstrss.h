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

#ifndef planstrss_h
#define planstrss_h

#include "nlstructuralelement.h"
#include "fei2dquadlin.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"

#include "directerrorindicatorrc.h"
#include "eleminterpmapperinterface.h"
#include "huertaerrorestimator.h"

namespace oofem {
/// Comment or uncomment the following line to force full or reduced integration
#define PlaneStress2d_reducedShearIntegration
/**
 * This class implements an isoparametric four-node quadrilateral plane-
 * stress elasticity finite element. Each node has 2 degrees of freedom.
 */
class PlaneStress2d : public NLStructuralElement, public ZZNodalRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
    public SpatialLocalizerInterface,
    public DirectErrorIndicatorRCInterface, public EIPrimaryUnknownMapperInterface,
    public HuertaErrorEstimatorInterface, public HuertaRemeshingCriteriaInterface
{
protected:
    static FEI2dQuadLin interpolation;
    int numberOfGaussPoints;

public:
    PlaneStress2d(int n, Domain *d);
    ~PlaneStress2d();

    virtual int computeNumberOfDofs(EquationID ut) { return 8; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);

    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }
    Interface *giveInterface(InterfaceType it);

    double computeVolumeAround(GaussPoint *gp);
    /**
     * Computes the global coordinates from given element's local coordinates.
     * Required by nonlocal material models.
     * @returns nonzero if successful
     */
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    /**
     * Computes the element local coordinates from given global coordinates.
     * @returns nonzero if successful (point inside); zero otherwise
     */
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    FEInterpolation *giveInterpolation() { return & interpolation; }

    Element *ZZNodalRecoveryMI_giveElement() { return this; }

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
    virtual double HuertaRemeshingCriteriaI_giveCharacteristicSize() { return DirectErrorIndicatorRCI_giveCharacteristicSize(); }
    virtual int HuertaRemeshingCriteriaI_givePolynOrder() { return 1; };

    virtual double DirectErrorIndicatorRCI_giveCharacteristicSize();

    virtual int EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType u,
                                                                 TimeStep *stepN, const FloatArray &coords,
                                                                 FloatArray &answer);
    virtual void EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer);

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void drawScalar(oofegGraphicContext &context);
    virtual void drawSpecial(oofegGraphicContext &);
    //void drawInternalState(oofegGraphicContext &);
#endif

    // definition & identification
    const char *giveClassName() const { return "PlaneStress2d"; }
    classType giveClassID() const { return PlaneStress2dClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_quad_1; }
    integrationDomain giveIntegrationDomain() { return _Square; }
    MaterialMode giveMaterialMode() { return _PlaneStress; }

protected:
    // edge load support
    void computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);

    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    void computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *gp, int i);
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);

    void computeGaussPoints();

    int giveApproxOrder() { return 1; }
    int giveNumberOfIPForMassMtrxIntegration() { return 4; }
};
} // end namespace oofem
#endif // planstrss_h
