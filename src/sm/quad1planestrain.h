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

#ifndef quad1planestrain_h
#define quad1planestrain_h

#include "structuralelement.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"
#include "directerrorindicatorrc.h"
#include "eleminterpmapperinterface.h"
#include "huertaerrorestimator.h"
#include "fei2dquadlin.h"

namespace oofem {
/// Comment or uncomment the following line to force full or reduced integration
///@todo Removed for now.
//#define Quad1PlaneStrain_reducedShearIntegration
/**
 * This class implements an isoparametric four-node quadrilateral plane-
 * stress structural finite element. Each node has 2 degrees of freedom.
 */
class Quad1PlaneStrain : public StructuralElement, public ZZNodalRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
    public SpatialLocalizerInterface,
    public DirectErrorIndicatorRCInterface, public EIPrimaryUnknownMapperInterface,
    public HuertaErrorEstimatorInterface, public HuertaRemeshingCriteriaInterface
{
protected:
    static FEI2dQuadLin interp;
    int numberOfGaussPoints;

public:
    Quad1PlaneStrain(int n, Domain *d);
    virtual ~Quad1PlaneStrain();

    virtual int computeNumberOfDofs(EquationID ut) { return 8; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    virtual FEInterpolation *giveInterpolation() { return &interp; }

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);

    virtual int testElementExtension(ElementExtension ext) { return ext == Element_EdgeLoadSupport; }

    virtual Interface *giveInterface(InterfaceType it);

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)  { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

    virtual double DirectErrorIndicatorRCI_giveCharacteristicSize();

    virtual int EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType u,
                                                                 TimeStep *stepN, const FloatArray &coords,
                                                                 FloatArray &answer);
    virtual void EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer);

    // HuertaErrorEstimatorInterface
    virtual Element *HuertaErrorEstimatorI_giveElement() { return this; }

    virtual void HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                                  IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                                  HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                                  int &localNodeId, int &localElemId, int &localBcId,
                                                                  IntArray &controlNode, IntArray &controlDof,
                                                                  HuertaErrorEstimator :: AnalysisMode aMode);
    virtual void HuertaErrorEstimatorI_computeLocalCoords(FloatArray &answer, const FloatArray &coords)
    { computeLocalCoordinates(answer, coords); }
    virtual void HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
    { computeNmatrixAt(aGaussPoint, answer); }

    // HuertaRemeshingCriteriaInterface
    virtual double HuertaRemeshingCriteriaI_giveCharacteristicSize() { return DirectErrorIndicatorRCI_giveCharacteristicSize(); };
    virtual int HuertaRemeshingCriteriaI_givePolynOrder() { return 1; };

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void drawScalar(oofegGraphicContext &context);
    virtual void drawSpecial(oofegGraphicContext &);
    //void drawInternalState(oofegGraphicContext &);
#endif

    // definition & identification
    virtual const char *giveClassName() const { return "Quad1PlaneStrain"; }
    virtual classType giveClassID() const { return Quad1PlaneStrainClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual Element_Geometry_Type giveGeometryType() const { return EGT_quad_1; }
    virtual integrationDomain giveIntegrationDomain() { return _Square; }
    virtual MaterialMode giveMaterialMode() { return _PlaneStrain; }

protected:
    // edge load support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *gp);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);

    virtual void computeGaussPoints();

    virtual int giveApproxOrder() { return 1; }
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 4; }
};
} // end namespace oofem
#endif // quad1planestrain_h
