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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef lspace_h
#define lspace_h

#include "sm/Elements/structural3delement.h"
#include "sm/ErrorEstimators/huertaerrorestimator.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "spatiallocalizer.h"

#define _IFT_LSpace_Name "lspace"

namespace oofem {
class FEI3dHexaLin;
class ParamKey;

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
class LSpace  : public Structural3DElement, public ZZNodalRecoveryModelInterface,
    public SPRNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface,
    public SpatialLocalizerInterface,
    public HuertaErrorEstimatorInterface
{
protected:
    static FEI3dHexaLin interpolation;
    bool reducedShearIntegration;
    static ParamKey IPK_LSpace_reducedShearIntegration; ///< [optional] Reduced shear integration flag.
public:
    LSpace(int n, Domain *d);
    virtual ~LSpace() { }
    FEInterpolation *giveInterpolation() const override;

    Interface *giveInterface(InterfaceType it) override;
    int testElementExtension(ElementExtension ext) override
    { return ( ( ( ext == Element_EdgeLoadSupport ) || ( ext == Element_SurfaceLoadSupport ) ) ? 1 : 0 ); }


    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override;
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override;

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep) override;

    void HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                          IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                          HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                          int &localNodeId, int &localElemId, int &localBcId,
                                                          IntArray &controlNode, IntArray &controlDof,
                                                          HuertaErrorEstimator :: AnalysisMode aMode) override;
    void HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_LSpace_Name; }
    const char *giveClassName() const override { return "LSpace"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_hexa_1;}

    void initializeFrom(InputRecord &ir, int priority) override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawSpecial(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawTriad(FloatArray &, int isurf);
#endif

protected:
    int giveNumberOfIPForMassMtrxIntegration() override { return 8; }
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui) override;
    /**
     * @name Surface load support
     */
    //@{
    int computeLoadLSToLRotationMatrix(FloatMatrix &answer, int iSurf, GaussPoint *gp) override;
    //@}
};
} // end namespace oofem
#endif // lspace_h
