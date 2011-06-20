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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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
     * @param n number associated with receiver
     * @param domain reference to domain.
     * @param startIndx first component, for which rule applies
     * @param endIndx last component, for which rule applies
     * @param dynamic flag indicating that receiver can change
     */
    LobattoIntegrationRule(int, Element *, int, int, bool dynamic);
    /// Destructor
    ~LobattoIntegrationRule();

    classType giveClassID() const { return LobattoIntegrationRuleClass; }
    const char *giveClassName() const { return "LobattoIntegrationRule"; }
    IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }

    int getRequiredNumberOfIntegrationPoints(integrationDomain dType, int approxOrder);

protected:
    int SetUpPointsOnLine(int, Element *, MaterialMode mode, GaussPoint ***gp);
    int SetUpPointsOnTriagle(int, Element *, MaterialMode mode, GaussPoint ***gp);
    int SetUpPointsOnSquare(int, Element *, MaterialMode mode, GaussPoint ***gp);
    int SetUpPointsOnCube(int, Element *, MaterialMode mode, GaussPoint ***gp);
    int SetUpPointsOnTetrahedra(int, Element *, MaterialMode mode, GaussPoint ***gp);
};
} // end namespace oofem
#endif // lobattoir_h
