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
#include "gausspoint.h"

#define _IFT_Quad1PlaneStrain_Name "quad1planestrain"

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

public:
    Quad1PlaneStrain(int n, Domain * d);
    virtual ~Quad1PlaneStrain();

    virtual int computeNumberOfDofs() { return 8; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    virtual FEInterpolation *giveInterpolation() const { return & interp; }

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);

    virtual int testElementExtension(ElementExtension ext) { return ext == Element_EdgeLoadSupport; }

    virtual Interface *giveInterface(InterfaceType it);

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

    virtual double DirectErrorIndicatorRCI_giveCharacteristicSize();

    virtual int EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType u,
                                                                 TimeStep *tStep, const FloatArray &coords,
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
    virtual void HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
    { computeNmatrixAt(* ( gp->giveLocalCoordinates() ), answer); }

    // HuertaRemeshingCriteriaInterface
    virtual double HuertaRemeshingCriteriaI_giveCharacteristicSize() { return DirectErrorIndicatorRCI_giveCharacteristicSize(); };
    virtual int HuertaRemeshingCriteriaI_givePolynOrder() { return 1; };

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &);
    virtual void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void drawScalar(oofegGraphicContext &context);
    virtual void drawSpecial(oofegGraphicContext &);
#endif

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_Quad1PlaneStrain_Name; }
    virtual const char *giveClassName() const { return "Quad1PlaneStrain"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _PlaneStrain; }

protected:
    // edge load support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *gp);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);

    virtual void computeGaussPoints();

    virtual int giveApproxOrder() { return 1; }
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 4; }
};
} // end namespace oofem
#endif // quad1planestrain_h
