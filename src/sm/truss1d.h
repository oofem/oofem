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

#ifndef truss1d_h
#define truss1d_h

#include "structuralelement.h"
#include "gaussintegrationrule.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "spatiallocalizer.h"
#include "directerrorindicatorrc.h"
#include "eleminterpmapperinterface.h"
#include "zzerrorestimator.h"
#include "mmashapefunctprojection.h"
#include "huertaerrorestimator.h"

namespace oofem {

/**
 * This class implements a two-node truss bar element for one-dimensional
 * analysis.
 *
 * A truss bar element is characterized by its 'length' and its 'pitch'. The
 * pitch is the angle in radians between the X-axis and the axis of the
 * element (oriented node1 to node2).
 * Note: element is formulated in global c.s.
 * Tasks:
 * - Calculating its Gauss points.
 * - Calculating its B,D,N matrices and dV.
 * - Expressing M,K,f,etc, in global axes. Methods like 'computeStiffness-
 *   Matrix' of class Element are here overloaded in order to account for
 *   rotational effects.
 */
class Truss1d : public StructuralElement,
    public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface, public SpatialLocalizerInterface,
    public DirectErrorIndicatorRCInterface,
    public EIPrimaryUnknownMapperInterface,
    public ZZErrorEstimatorInterface, public ZZRemeshingCriteriaInterface, public MMAShapeFunctProjectionInterface,
    public HuertaErrorEstimatorInterface, public HuertaRemeshingCriteriaInterface
{
protected:
    double length;

public:
    Truss1d(int n, Domain *d);
    ~Truss1d() { }

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    virtual int computeNumberOfDofs(EquationID ut) { return 2; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    // characteristic length in gp (for some material models)
    double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
    { return this->giveLength(); }

    double computeVolumeAround(GaussPoint *gp);

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    Interface *giveInterface(InterfaceType it);

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    void drawScalar(oofegGraphicContext &context);
#endif

    // definition & identification
    const char *giveClassName() const { return "Truss1d"; }
    classType giveClassID() const { return Truss1dClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

    // ZZNodalRecoveryMInterface
    //void ZZNodalRecoveryMI_computeNValProduct (FloatArray& answer, InternalStateType type, TimeStep* tStep);
    //void ZZNodalRecoveryMI_computeNNMatrix (FloatArray& answer, InternalStateType type);
    int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    Element *ZZNodalRecoveryMI_giveElement() { return this; }
    void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint,
                                                             InternalStateType type);

    // NodalAveragingRecoveryMInterface
    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);

    void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }

    // SpatialLocalizerInterface
    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

    virtual double DirectErrorIndicatorRCI_giveCharacteristicSize();

    // EIPrimaryUnknownMInterface
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
    virtual double HuertaRemeshingCriteriaI_giveCharacteristicSize() { return DirectErrorIndicatorRCI_giveCharacteristicSize(); }
    virtual int HuertaRemeshingCriteriaI_givePolynOrder() { return 1; };

    virtual void MMAShapeFunctProjectionInterface_interpolateIntVarAt(FloatArray &answer, FloatArray &coords,
                                                                      coordType ct, nodalValContainerType &list,
                                                                      InternalStateType type, TimeStep *tStep);

    integrationDomain giveIntegrationDomain() { return _Line; }
    MaterialMode giveMaterialMode() { return _1dMat; }

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeGaussPoints();

    double giveLength();
    double givePitch();
    int giveApproxOrder() { return 1; }
};
} // end namespace oofem
#endif // truss1d_h
