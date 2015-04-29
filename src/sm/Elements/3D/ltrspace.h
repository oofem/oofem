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

#ifndef ltrspace_h
#define ltrspace_h

#include "Elements/structural3delement.h"
#include "ErrorEstimators/directerrorindicatorrc.h"
#include "ErrorEstimators/huertaerrorestimator.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"
#include "ErrorEstimators/zzerrorestimator.h"
#include "mmashapefunctprojection.h"

#define _IFT_LTRSpace_Name "ltrspace"

namespace oofem {
class FEI3dTetLin;

/**
 * This class implements a linear tetrahedral four-node finite element for stress analysis.
 * Each node has 3 degrees of freedom.
 */
class LTRSpace : public Structural3DElement, public ZZNodalRecoveryModelInterface,
public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
public SpatialLocalizerInterface,
public ZZErrorEstimatorInterface,
public HuertaErrorEstimatorInterface
{
protected:
    static FEI3dTetLin interpolation;

public:
    LTRSpace(int n, Domain * d);
    virtual ~LTRSpace() { }

    virtual FEInterpolation *giveInterpolation() const;

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 4; }
    virtual Interface *giveInterface(InterfaceType it);

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_LTRSpace_Name; }
    virtual const char *giveClassName() const { return "LTRSpace"; }
    
#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawSpecial(oofegGraphicContext &gc, TimeStep *tStep);
#endif

    virtual double giveRelativeSelfComputationalCost() { return 2.15; }

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep);

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

    // HuertaErrorEstimatorInterface
    virtual void HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                                  IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                                  HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                                  int &localNodeId, int &localElemId, int &localBcId,
                                                                  IntArray &controlNode, IntArray &controlDof,
                                                                  HuertaErrorEstimator :: AnalysisMode aMode);
    virtual void HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);

protected:

    /**
     * @name Surface load support
     */
    //@{
    virtual IntegrationRule *GetSurfaceIntegrationRule(int);
    //@}
};
} // end namespace oofem
#endif // ltrspace_h
