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

#ifndef patchintegrationrule_h
#define patchintegrationrule_h

#include "gaussintegrationrule.h"

namespace oofem {
class FEI2dTrLin;
class FEI3dTrQuad;
class Triangle;

/**
 * PatchIntegrationRule provides integration over a triangle patch.
 * Input to the constructor is:
 *  -int n:         number of quadrature points per triangle.
 *  -Element *e:    parent element pointer
 *  -iTriangles:    array of triangles describing the subdivision of the element.
 *
 * @author Erik Svenning (Major modifications)
 *
 */
class OOFEM_EXPORT PatchIntegrationRule : public GaussIntegrationRule
{
protected:
    std :: vector< Triangle >mTriangles;

    // Interpolation used to distribute quadrature points
    // in each triangle of the patch.
    static FEI2dTrLin mTriInterp;
    static FEI3dTrQuad mTriInterpQuad;


public:
    /// Constructor.
    PatchIntegrationRule(int n, Element *e, const std :: vector< Triangle > &iTriangles);
    /// Destructor.
    virtual ~PatchIntegrationRule();

    virtual const char *giveClassName() const { return "PatchIntegrationRule"; }

    // TODO: Give this function a better name.
    // Note: the fact that this function is inherited complicates name change.
    virtual int SetUpPointsOnTriangle(int nPoints, MaterialMode mode);
    virtual int SetUpPointsOnWedge(int nPointsTri, int nPointsDepth, MaterialMode mode);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj);
};
} // end namespace oofem
#endif // patchintegrationrule_h
