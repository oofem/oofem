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

#ifndef tria2platesubsoil_H
#define tria2platesubsoil_H

#include "sm/Elements/tria1platesubsoil.h"

#define _IFT_Tria2PlateSubSoil_Name "tria2platesubsoil"

namespace oofem {
class FEI2dTrQuad;

/**
 * This class implements an triangular six-node plate subsoil element with quadratic interpolation in xy plane.
 * Each node has 1 degree of freedom (out-of-plane displacement).
 *
 * Loading types supported;
 *
 * Reference:
 * Bittnar, Sejnoha: Numerical Methods in Structural Mechanics, Thomas Telford, Jan 1, 1996, ISBN-13: 978-0784401705
 * @author Borek Patzak
 */
class Tria2PlateSubSoil : public Tria1PlateSubSoil
{
protected:
    static FEI2dTrQuad interp_quad;

public:
    Tria2PlateSubSoil(int n, Domain * d);
    virtual ~Tria2PlateSubSoil() { }

    FEInterpolation *giveInterpolation() const override;
    FEInterpolation *giveInterpolation(DofIDItem id) const override;
    
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_Tria2PlateSubSoil_Name; }
    const char *giveClassName() const override { return "Tria2PlateSubSoil"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_triangle_2;}

    int computeNumberOfDofs() override { return 6; }
    
    
protected:
    void computeGaussPoints() override;
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override;
    
    /**
     * @name Surface load support
     */
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;    
    void computeSurfaceNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords) override;
};
} // end namespace oofem
#endif // tria2platesubsoil_H
