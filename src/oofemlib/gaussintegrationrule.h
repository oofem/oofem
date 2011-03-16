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

#ifndef gaussintegrationrule_h
#define gaussintegrationrule_h

#include "integrationrule.h"
#include "element.h"

namespace oofem {
/**
 * Class representing Gaussian-quadrature integration rule.
 * The number of integration points and their coordinates and integration weights depends on
 * integration rule type (rule for integration in 1d, 2d, 3d) and required  acurracy.
 */
class GaussIntegrationRule : public IntegrationRule
{
    /*
     * DESCRIPTION:
     * Implements integration rule class.
     * Stores integration points used for integration
     * of necesary terms (for example computation of  stiffness matrix
     * or computation of element nodal force vector )
     * and it  corresponds to some local strains
     * on finite element level. Finite element can have many
     * integration rules corresponding to  different strains.
     *
     * TASKS:
     * instanciating yourself
     * returning number of integration points used
     * returning requested integration point - method getIntegrationPoint
     * returning inteval of components (i.e.of local strain vector), where apply
     * printing yourself
     * updating yourself
     * initializing for new time step
     * saving & restoring context
     */
public:
    /**
     * Constructor.
     * @param n number associated with receiver
     * @param domain reference to domain.
     * @param startIndx first component, for which rule applies
     * @param endIndx last component, for which rule applies
     * @param dynamic flag indicating that receiver can change
     */
    GaussIntegrationRule(int n, Element *e, int startIndx, int endIndx, bool dynamic = false);
    GaussIntegrationRule(int n, Element *e);
    /// Destructor
    virtual ~GaussIntegrationRule();

    ///Returns classType id of receiver.
    classType giveClassID() const { return GaussIntegrationRuleClass; }
    ///Returns class name of the receiver.
    const char *giveClassName() const { return "GaussIntegrationRule"; }
    IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }

    /**
     * Returns requred number of integration points to exactly integrate
     * polynomial of order approxOrder on given domain.
     * When approxOrder is too large and is not supported by implementation
     * method returns -1.
     */
    int getRequiredNumberOfIntegrationPoints(integrationDomain dType, int approxOrder);

protected:
    /**
     * Sets up receiver's  integration points on unit line integration domain.
     * @returns number of integration points.
     */
    int SetUpPointsOnLine(int, MaterialMode, GaussPoint * * *);
    /**
     * Sets up receiver's  integration points on triangular (area coords) integration domain.
     * @returns number of integration points.
     */
    virtual int SetUpPointsOnTriagle(int, MaterialMode, GaussPoint * * *);
    /**
     * Sets up receiver's  integration points on unit square integration domain.
     * @returns number of integration points.
     */
    int SetUpPointsOnSquare(int, MaterialMode, GaussPoint * * *);
    /**
     * Sets up receiver's  integration points on unit cube integration domain.
     * @returns number of integration points.
     */
    int SetUpPointsOnCube(int, MaterialMode, GaussPoint * * *);
    /**
     * Sets up receiver's  integration points on tetrahedra (volume coords) integration domain.
     * @returns number of integration points.
     */
    int SetUpPointsOnTetrahedra(int, MaterialMode, GaussPoint * * *);
    /**
     * Sets up integration points on 2D embedded line inside 2D volume (the list of local coordinates
     * should be provided).
     */
    int SetUpPointsOn2DEmbeddedLine(int nPoints, MaterialMode mode, GaussPoint ***,
                                    const FloatArray **coords);
};
} // end namespace oofem
#endif // gaussintegrationrule_h
