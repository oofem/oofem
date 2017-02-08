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

#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/quad1platesubsoil.h"

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

    virtual FEInterpolation *giveInterpolation() const;
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const;

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_Quad2PlateSubSoil_Name; }
    virtual const char *giveClassName() const { return "Quad2PlateSubSoil"; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int computeNumberOfDofs() { return 8; }

protected:
    virtual void computeGaussPoints();
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
};
} // end namespace oofem
#endif // quad2platesubsoil_H
