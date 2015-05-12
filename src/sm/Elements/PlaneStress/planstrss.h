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

#ifndef planstrss_h
#define planstrss_h

#include "Elements/structural2delement.h"
#include "ErrorEstimators/directerrorindicatorrc.h"
#include "ErrorEstimators/huertaerrorestimator.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"

#define _IFT_PlaneStress2d_Name "planestress2d"

namespace oofem {
class FEI2dQuadLin;

/// Comment or uncomment the following line to force full or reduced integration
#define PlaneStress2d_reducedShearIntegration
/**
 * This class implements an isoparametric four-node quadrilateral plane-
 * stress elasticity finite element. Each node has 2 degrees of freedom.
 */
class PlaneStress2d : public PlaneStressElement, public ZZNodalRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
public SpatialLocalizerInterface,
public HuertaErrorEstimatorInterface
{
protected:
    static FEI2dQuadLin interpolation;

public:
    PlaneStress2d(int n, Domain * d);
    virtual ~PlaneStress2d();

    virtual double giveCharacteristicSize(GaussPoint *gp, FloatArray &normalToCrackPlane, ElementCharSizeMethod method);
    virtual double giveParentElSize() const { return 4.0; }

    virtual Interface *giveInterface(InterfaceType it);
    virtual FEInterpolation *giveInterpolation() const;

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual void HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                                  IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                                  HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                                  int &localNodeId, int &localElemId, int &localBcId,
                                                                  IntArray &controlNode, IntArray &controlDof,
                                                                  HuertaErrorEstimator :: AnalysisMode aMode);
    virtual void HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawSpecial(oofegGraphicContext &gc, TimeStep *tStep);
#endif

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_PlaneStress2d_Name; }
    virtual const char *giveClassName() const { return "PlaneStress2d"; }
    virtual IRResultType initializeFrom(InputRecord *ir);

protected:

    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);

    int giveNumberOfIPForMassMtrxIntegration() { return 4; } 
};
} // end namespace oofem
#endif // planstrss_h
