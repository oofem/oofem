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

#ifndef fei2dquadbiquad_h
#define fei2dquadbiquad_h

#include "fei2dquadquad.h"

namespace oofem {
/**
 * Class representing a 2d quadrilateral with bi-quadratic interpolation based on isoparametric coordinates.
 * Local Node Numbering
 *       ^ eta
 *       |
 * (4)--(7)--(3)
 *  |         |
 *  |         |
 * (8)  (9)  (6)-->ksi
 *  |         |
 *  |         |
 * (1)--(5)--(2)
 * Everything regarding edges can be directly inherited by FEI2dQuadQuad.
 * @note Untested.
 * @author Mikael Öhman
 */
class OOFEM_EXPORT FEI2dQuadBiQuad : public FEI2dQuadQuad
{
public:
    FEI2dQuadBiQuad(int ind1, int ind2) : FEI2dQuadQuad(ind1, ind2) { }

    // Bulk
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int giveNumberOfNodes() const { return 9; }
    
    virtual IntegrationRule *giveIntegrationRule(int order);

protected:
    virtual void giveDerivatives(FloatMatrix &answer, const FloatArray &lcoords);
};
} // end namespace oofem
#endif // fei2dquadbiquad_h
