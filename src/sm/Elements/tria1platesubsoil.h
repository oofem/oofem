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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#ifndef tria1platesubsoil_H
#define tria1platesubsoil_H

#include "sm/Elements/structuralelement.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

#define _IFT_Tria1PlateSubSoil_Name "tria1platesubsoil"

namespace oofem {
class FEI2dTrLin;

/**
 * This class implements an triangular four-node plate subsoil element with linear interpolation in xy plane.
 * Each node has 1 degree of freedom (out-of-plane displacement).
 *
 * Loading types supported;
 *
 * Reference:
 * Bittnar, Sejnoha: Numerical Methods in Structural Mechanics, Thomas Telford, Jan 1, 1996, ISBN-13: 978-0784401705
 * @author Borek Patzak
 */
class Tria1PlateSubSoil : public StructuralElement,
public ZZNodalRecoveryModelInterface,
public SPRNodalRecoveryModelInterface
{
protected:
    static FEI2dTrLin interp_lin;

public:
    Tria1PlateSubSoil(int n, Domain * d);
    virtual ~Tria1PlateSubSoil() { }

    FEInterpolation *giveInterpolation() const override;
    FEInterpolation *giveInterpolation(DofIDItem id) const override;

    MaterialMode giveMaterialMode() override { return _2dPlateSubSoil; }
    int testElementExtension(ElementExtension ext) override { return ( ( ext == Element_SurfaceLoadSupport ) ? 1 : 0 ); }

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_Tria1PlateSubSoil_Name; }
    const char *giveClassName() const override { return "Tria1PlateSubSoil"; }
    void initializeFrom(InputRecord &ir) override;

    int computeNumberOfDofs() override { return 3; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;
    Element_Geometry_Type giveGeometryType() const override {return EGT_triangle_1;}


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

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override { return this->numberOfGaussPoints; }
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override { return SPRPatchType_2dxy; }
    /**
     * @name Surface load support
     */
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;    
    void computeSurfaceNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords) override;
    
    //@{
    //void computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp);
    //void giveSurfaceDofMapping(IntArray &answer, int iSurf) const override;
    //double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf) override;
    //int computeLoadLSToLRotationMatrix(FloatMatrix &answer, int iSurf, GaussPoint *gp) override;
    //@}
};
} // end namespace oofem
#endif // tria1platesubsoil_H
