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

#ifndef truss3d_h
#define truss3d_h

#include "sm/Elements/nlstructuralelement.h"
#include "sm/ErrorEstimators/directerrorindicatorrc.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"

#define _IFT_Truss3d_Name "truss3d"

namespace oofem {
class FEI3dLineLin;

/**
 * This class implements a two-node truss bar element for three-dimensional
 * analysis.
 */
class Truss3d : public NLStructuralElement,
    public ZZNodalRecoveryModelInterface,
    public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI3dLineLin interp;

public:
    Truss3d(int n, Domain *d);
    virtual ~Truss3d() { }

    FEInterpolation *giveInterpolation() const override;

    double computeLength() override;

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { this->computeLumpedMassMatrix(answer, tStep); }
    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    int computeNumberOfDofs() override { return 6; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;


    // characteristic length (for crack band approach)
    double giveCharacteristicLength(const FloatArray &normalToCrackPlane) override
    { return this->computeLength(); }

    double computeVolumeAround(GaussPoint *gp) override;

    int testElementExtension(ElementExtension ext) override { return ( ext == Element_EdgeLoadSupport ); }

    Interface *giveInterface(InterfaceType it) override;

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep) override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
#endif

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_Truss3d_Name; }
    const char *giveClassName() const override { return "Truss3d"; }
    void initializeFrom(InputRecord &ir) override;
    MaterialMode giveMaterialMode() override { return _1dMat; }
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

protected:
    // edge load support
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const override;
    double computeEdgeVolumeAround(GaussPoint *gp, int) override;
    int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *gp) override;
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;
    void computeGaussPoints() override;
};
} // end namespace oofem
#endif // truss3d_h
