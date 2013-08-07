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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef patchintegrationrule_h
#define patchintegrationrule_h

#define PATCH_INT_DEBUG 1

#include "gaussintegrationrule.h"
#include "classtype.h"
//#include "patch.h"

namespace oofem {
class FEI2dTrLin;
class Triangle;

/**
 * PatchIntegrationRule provides integration over a triangle patch.
 * Input to the constructor is:
 * 	-int n: 		number of quadrature points per triangle.
 * 	-Element *e:	parent element pointer
 * 	-iTriangles:	array of triangles describing the subdivision of the element.
 *
 *  * @author Erik Svenning (Major modifications)
 *
 */
class PatchIntegrationRule : public GaussIntegrationRule
{
protected:
	std::vector<Triangle> mTriangles;

	// Interpolation used to distribute quadrature points
	// in each triangle of the patch.
    static FEI2dTrLin mTriInterp;


public:
    /// Constructor.
    PatchIntegrationRule(int n, Element *e, const std::vector<Triangle> &iTriangles);
    /// Destructor.
    virtual ~PatchIntegrationRule();

    // TODO: Give this function a better name.
    // Note: the fact that this function is inherited complicates name change.
    virtual int SetUpPointsOnTriangle(int nPoints, MaterialMode mode);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj);
    virtual classType giveClassID() const { return PatchIntegrationRuleClass; }
};
} // end namespace oofem
#endif // patchintegrationrule_h
