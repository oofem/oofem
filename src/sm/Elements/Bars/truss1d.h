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

#include "sm/Elements/nlstructuralelement.h"
#include "sm/ErrorEstimators/directerrorindicatorrc.h"
#include "sm/ErrorEstimators/zzerrorestimator.h"
#include "sm/ErrorEstimators/huertaerrorestimator.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "spatiallocalizer.h"
#include "mmashapefunctprojection.h"

#define _IFT_Truss1d_Name "truss1d"

namespace oofem {
class FEI1dLin;

/**
 * This class implements a two-node truss bar element for one-dimensional
 * analysis.
 */
class Truss1d : public NLStructuralElement,
    public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface, public SpatialLocalizerInterface,
    public ZZErrorEstimatorInterface,
    public HuertaErrorEstimatorInterface
{
protected:
    static FEI1dLin interp;

public:
    Truss1d(int n, Domain *d);
    virtual ~Truss1d() { }

    FEInterpolation *giveInterpolation() const override;

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { computeLumpedMassMatrix(answer, tStep); }

    int computeNumberOfDofs() override { return 2; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;

    // characteristic length (for crack band approach)
    double giveCharacteristicLength(const FloatArray &normalToCrackPlane) override
    { return this->computeLength(); }

    double computeVolumeAround(GaussPoint *gp) override;

    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

    int testElementExtension(ElementExtension ext) override { return 0; }

    Interface *giveInterface(InterfaceType it) override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_Truss1d_Name; }
    const char *giveClassName() const override { return "Truss1d"; }
    MaterialMode giveMaterialMode() override { return _1dMat; }

    // NodalAveragingRecoveryMInterface
    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep) override;

    // HuertaErrorEstimatorInterface
    void HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                          IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                          HuertaErrorEstimatorInterface::SetupMode sMode, TimeStep *tStep,
                                                          int &localNodeId, int &localElemId, int &localBcId,
                                                          IntArray &controlNode, IntArray &controlDof,
                                                          HuertaErrorEstimator::AnalysisMode aMode) override;
    void HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override; ///@todo BH is just used in nlstructuralelement
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;
    void computeGaussPoints() override;
};
} // end namespace oofem
#endif // truss1d_h
