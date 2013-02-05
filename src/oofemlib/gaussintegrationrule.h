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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef gaussintegrationrule_h
#define gaussintegrationrule_h

#include "integrationrule.h"
#include "element.h"

namespace oofem {
/**
 * Class representing Gaussian-quadrature integration rule.
 * The number of integration points and their coordinates and integration weights depends on
 * integration rule type (rule for integration in 1d, 2d, 3d) and required accuracy.
 * The positions and weights are determined by the minimum required of points to integrate a polynomial exactly (while the points are strictly within the domain)
 * 
 * Tasks:
 * - Returning number of integration points used
 * - Returning requested integration point
 * - Updating itself
 * - Saving and restoring context
 * 
 * @see GaussPoint
 */
class GaussIntegrationRule : public IntegrationRule
{
public:
    /**
     * Constructor.
     * @param n Number associated with receiver.
     * @param e Element associated with receiver.
     * @param startIndx First component, for which rule applies.
     * @param endIndx Last component, for which rule applies.
     * @param dynamic Flag indicating that receiver can change.
     */
    GaussIntegrationRule(int n, Element *e, int startIndx, int endIndx, bool dynamic = false);
    GaussIntegrationRule(int n, Element *e);
    /// Destructor
    virtual ~GaussIntegrationRule();

    virtual classType giveClassID() const { return GaussIntegrationRuleClass; }
    virtual const char *giveClassName() const { return "GaussIntegrationRule"; }
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }

    virtual int getRequiredNumberOfIntegrationPoints(integrationDomain dType, int approxOrder);


    virtual int SetUpPointsOnLine(int, MaterialMode);
    virtual int SetUpPointsOnTriangle(int, MaterialMode);
    virtual int SetUpPointsOnSquare(int, MaterialMode);
    virtual int SetUpPointsOnCube(int, MaterialMode);
    virtual int SetUpPointsOnTetrahedra(int, MaterialMode);

    virtual int SetUpPointsOn2DEmbeddedLine(int nPoints, MaterialMode mode, const FloatArray **coords);

    void giveTriCoordsAndWeights(int nPoints, FloatArray &coords_xi1, FloatArray &coords_xi2, FloatArray &weights);
    void giveLineCoordsAndWeights(int nPoints, FloatArray &coords_xi, FloatArray &weights);

    virtual int SetUpPointsOnWedge(int nPointsTri, int nPointsDepth, MaterialMode mode);

};
} // end namespace oofem
#endif // gaussintegrationrule_h
