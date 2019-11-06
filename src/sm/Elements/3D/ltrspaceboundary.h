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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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

#ifndef ltrspaceboundary_h
#define ltrspaceboundary_h

#include "sm/Elements/structural3delement.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"
#include "mmashapefunctprojection.h"

#define _IFT_LTRSpaceBoundary_Name "ltrspaceboundary"
#define _IFT_LTRSpaceBoundary_Location "location"

namespace oofem {
class FEI3dTetLin;

/**
 * This class implements a linear tetrahedral four-node finite element.
 * Each node has 3 degrees of freedom. This element is used for 3D RVE analyses with Periodic Boundary Conditions.
 * At least one node is located at the image boundary.
 * These nodes are replaced with a periodic mirror nodes and a control node is used to impose the macroscopic (average) strain.
 * MACROSCOPIC INPUT: DEFORMATION GRADIENT TENSOR (3D, 9 COMPONENTS: Exx Exy Exz Eyx Eyy Eyz Ezx Ezy Ezz)
 *
 * @author: Adam Sciegaj
 */
class LTRSpaceBoundary : public Structural3DElement, public NodalAveragingRecoveryModelInterface,
    public SpatialLocalizerInterface
{
protected:
    static FEI3dTetLin interpolation;
    IntArray location;
    virtual void computeTransformationMatrix(FloatMatrix &answer, TimeStep *tStep);
    void giveSwitches(IntArray &answer, int location);

    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override;
    void computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    double computeVolumeAround(GaussPoint *gp) override;
    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    double giveLengthInDir(const FloatArray &normalToCrackPlane) override;

public:
    LTRSpaceBoundary(int n, Domain *d);
    virtual ~LTRSpaceBoundary() { }

    FEInterpolation *giveInterpolation() const override;

    Interface *giveInterface(InterfaceType it) override;

    int computeNumberOfDofs() override { return 21; };
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    const IntArray giveLocation() override { return location; };
    void recalculateCoordinates(int nodeNumber, FloatArray &coords) override;

    // definition & identification
    void initializeFrom(InputRecord &ir) override;
    const char *giveInputRecordName() const override { return _IFT_LTRSpaceBoundary_Name; }
    const char *giveClassName() const override { return "LTRSpaceBoundary"; }

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // LTRSpaceBoundary_h
