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
#include "fei1dlin.h"
#include "gausspoint.h"

#define _IFT_Truss1d_Name "truss1d"

namespace oofem {
/**
 * This class implements a two-node truss bar element for one-dimensional
 * analysis.
 */
class Truss1d : public StructuralElement,
public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface, public SpatialLocalizerInterface,
public EIPrimaryUnknownMapperInterface,
public ZZErrorEstimatorInterface, public MMAShapeFunctProjectionInterface,
public HuertaErrorEstimatorInterface
{
protected:
    static FEI1dLin interp;

public:
    Truss1d(int n, Domain * d);
    virtual ~Truss1d() { }

    virtual FEInterpolation *giveInterpolation() const { return & interp; }

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }

    virtual int computeNumberOfDofs() { return 2; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    // characteristic length in gp (for some material models)
    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
    { return this->computeLength(); }

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    virtual Interface *giveInterface(InterfaceType it);

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &);
    virtual void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void drawScalar(oofegGraphicContext &context);
#endif

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_Truss1d_Name; }
    virtual const char *giveClassName() const { return "Truss1d"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _1dMat; }

    // ZZNodalRecoveryMInterface
    //void ZZNodalRecoveryMI_computeNValProduct (FloatArray& answer, InternalStateType type, TimeStep* tStep);
    //void ZZNodalRecoveryMI_computeNNMatrix (FloatArray& answer, InternalStateType type);

    // NodalAveragingRecoveryMInterface
    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep);

    // SpatialLocalizerInterface
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

    // EIPrimaryUnknownMInterface
    virtual void EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
                                                                       TimeStep *tStep, const FloatArray &lcoords,
                                                                       FloatArray &answer);
    virtual void EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer);

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
    { computeNmatrixAt(* ( gp->giveLocalCoordinates() ), answer); }

    virtual void MMAShapeFunctProjectionInterface_interpolateIntVarAt(FloatArray &answer, FloatArray &coords,
                                                                      coordType ct, nodalValContainerType &list,
                                                                      InternalStateType type, TimeStep *tStep);

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
    virtual void computeGaussPoints();
    virtual int giveApproxOrder() { return 1; }
};
} // end namespace oofem
#endif // truss1d_h
