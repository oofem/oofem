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

#ifndef fei3dhexatriquad_h
#define fei3dhexatriquad_h

#include "fei3dhexaquad.h"

namespace oofem {
/**
 * Class representing implementation of tri-quadratic hexahedra interpolation class.
 * @author Mikael Öhman
 */
class OOFEM_EXPORT FEI3dHexaTriQuad : public FEI3dHexaQuad
{
public:
    FEI3dHexaTriQuad() : FEI3dHexaQuad() { }

    virtual integrationDomain giveIntegrationDomain() const { return _Cube; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_hexa_27; }

    // Bulk
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int giveNumberOfNodes() const { return 27; }
    
    // Surface
    virtual void surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    //virtual void surfaceEvaldNdx(FloatMatrix&answer, int isurf, const FloatArray& lcoords, const FEICellGeometry& cellgeo);
    virtual void computeLocalSurfaceMapping(IntArray &nodes, int iSurf);
    virtual double evalNXIntegral(int iSurf, const FEICellGeometry &cellgeo);

    virtual IntegrationRule *giveIntegrationRule(int order);
    virtual IntegrationRule *giveBoundaryIntegrationRule(int order, int boundary);

protected:
    virtual void giveLocalDerivative(FloatMatrix &dN, const FloatArray &lcoords);
};
} // end namespace oofem
#endif
