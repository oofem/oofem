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
 *               Copyrighlint (C) 1993 - 2014   Borek Patzak
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

#ifndef quad2platesubsoil_H
#define quad2platesubsoil_H

#include "sm/Elements/structuralelement.h"
#include "sm/Elements/quad1platesubsoil.h"

#define _IFT_Quad2PlateSubSoil_Name "quad2platesubsoil"

namespace oofem {
class FEI2dQuadQuad;

/**
 * This class implements a quadrilateral eight-node plate subsoil element in xy plane.
 * Each node has 1 degree of freedom (out-of-plane displacement).
 *
 * Loading types supported;
 *
 * Reference:
 * Bittnar, Sejnoha: Numerical Methods in Structural Mechanics, Thomas Telford, Jan 1, 1996, ISBN-13: 978-0784401705
 * @author Borek Patzak
 */
class Quad2PlateSubSoil : public Quad1PlateSubSoil
{
protected:
    static FEI2dQuadQuad interp_quad;

public:
    Quad2PlateSubSoil(int n, Domain * d);
    virtual ~Quad2PlateSubSoil() { }

    FEInterpolation *giveInterpolation() const override;
    FEInterpolation *giveInterpolation(DofIDItem id) const override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_Quad2PlateSubSoil_Name; }
    const char *giveClassName() const override { return "Quad2PlateSubSoil"; }
    void initializeFrom(InputRecord &ir) override;

    int computeNumberOfDofs() override { return 8; }

protected:
    void computeGaussPoints() override;
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    /**
     * @name Surface load support
     */
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;    
    void computeSurfaceNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords) override;
    
};
} // end namespace oofem
#endif // quad2platesubsoil_H
