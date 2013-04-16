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

#ifndef lobattoir_h
#define lobattoir_h

#include "integrationrule.h"

namespace oofem {
/**
 * Class representing Lobatto-quadrature integration rule.
 * The number of integration points and their coordinates and integration weights depends on
 * integration rule type (rule for integration in 1d, 2d, 3d) and required  accuracy.
 */
class LobattoIntegrationRule : public IntegrationRule
{
public:
    /**
     * Constructor.
     * @param n Number associated with receiver
     * @param e Reference to engineering model.
     * @param startIndx First component, for which rule applies.
     * @param endIndx Last component, for which rule applies.
     * @param dynamic Flag indicating that receiver can change.
     */
    LobattoIntegrationRule(int n, Element *e, int startIndx, int endIndx, bool dynamic);
    /// Destructor
    virtual ~LobattoIntegrationRule();

    virtual classType giveClassID() const { return LobattoIntegrationRuleClass; }
    virtual const char *giveClassName() const { return "LobattoIntegrationRule"; }
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }

    virtual int getRequiredNumberOfIntegrationPoints(integrationDomain dType, int approxOrder);

    virtual int SetUpPointsOnLine(int, MaterialMode mode);
    virtual int SetUpPointsOnTriangle(int, MaterialMode mode);
    virtual int SetUpPointsOnSquare(int, MaterialMode mode);
    virtual int SetUpPointsOnCube(int, MaterialMode mode);
    virtual int SetUpPointsOnTetrahedra(int, MaterialMode mode);
};
} // end namespace oofem
#endif // lobattoir_h
