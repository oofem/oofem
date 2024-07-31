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

#ifndef quad1mindlin_H
#define quad1mindlin_H

#include "sm/Elements/structuralelement.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

#define _IFT_Quad1Mindlin_Name "quad1mindlin"
#define _IFT_Quad1Mindlin_ReducedIntegration "reducedintegration"

namespace oofem {
class FEI2dQuadLin;

/**
 * This class implements an quadrilateral four-node Mindlin plate.
 * Each node has 3 degrees of freedom (out-of-plane displacement, in-plane rotations).
 * This type of element exhibit strong shear locking (thin plates exhibit almost no bending).
 * No reduced integration is used, as it causes numerical problems.
 *
 * Loading types supported;
 * - Gravity load
 *
 * Reference:
 * Robert Cook, David Malkus, Michael Plesha
 * Concepts and Applications of Finite Element Analysis - Third edition
 * ISBN: 0-471-84788-7
 *
 * @author Mikael Öhman
 */
class Quad1Mindlin : public StructuralElement,
    public ZZNodalRecoveryModelInterface,
    public SPRNodalRecoveryModelInterface
{
protected:
    static FEI2dQuadLin interp_lin;
    /// Flag controlling reduced (one - point) integration for shear
    bool reducedIntegrationFlag = false;

public:
    Quad1Mindlin(int n, Domain *d);

    FEInterpolation *giveInterpolation() const override;
    FEInterpolation *giveInterpolation(DofIDItem id) const override;

    MaterialMode giveMaterialMode() override { return _2dPlate; }
    int testElementExtension(ElementExtension ext) override { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_Quad1Mindlin_Name; }
    const char *giveClassName() const override { return "Quad1Mindlin"; }
    void initializeFrom(InputRecord &ir) override;

    int computeNumberOfDofs() override { return 12; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;
    Element_Geometry_Type giveGeometryType() const override {return EGT_quad_1;}


    void computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp) override;

    double giveCharacteristicLength(const FloatArray &normalToCrackPlane) override;
    double computeVolumeAround(GaussPoint *gp) override;

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { computeLumpedMassMatrix(answer, tStep); }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    Interface *giveInterface(InterfaceType it) override;

protected:
    void computeGaussPoints() override;
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode) override;
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;

    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

    void giveEdgeDofMapping(IntArray &answer, int iEdge) const override;
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) override;
    int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp) override;

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override { return this->numberOfGaussPoints; }
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override { return SPRPatchType_2dxy; }
};
} // end namespace oofem
#endif // quad1mindlin_H
