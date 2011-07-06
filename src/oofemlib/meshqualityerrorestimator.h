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

#ifndef meshqualityerrorestimator_h
#define meshqualityerrorestimator_h

#include "errorestimator.h"
#include "classtype.h"
#include "errorestimatortype.h"

namespace oofem {
class Domain;
class Element;
class TimeStep;
class IntegrationRule;
class FEInterpolation;
/**
 * This error estimator measures the quality of the elements.
 * This is useful when elements are deformed during the simulation.
 *
 * Error estimates should be unitless, positive, and scaled properly.
 * Zero error should reflect a perfect element.
 * Error value of one should reflect an element of acceptable quality.
 *
 * @note Its open for discussion whether the routine should give the reciprocal error to handle extremely badly shaped elements.
 *
 * @author Mikael Öhman
 */
class MeshQualityErrorEstimator : public ErrorEstimator
{
protected:
    /**
     * Computes error based on the inscribed triangle/circle ratio.
     * @return Error value for the element. Zero for equilateral triangle.
     */
    static double computeTriangleRadiusError(Element *elem);

    /**
     * Computes the error based on the conditioning of the Jacobian.
     * @return Error value for the element.
     */
    static double computeJacobianError(FEInterpolation &fei, IntegrationRule &ir, Element *elem);

public:
    /// Constructor
    MeshQualityErrorEstimator(int n, Domain *d) : ErrorEstimator(n, d) { this->eeType = EET_MeshQuality; }
    /// Destructor
    ~MeshQualityErrorEstimator() { }

    virtual double giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep);

    /// Gives the max error from any element in the domain.
    virtual double giveValue(EE_ValueType type, TimeStep *tStep);

    /// Empty implementation.
    virtual int estimateError(EE_ErrorMode mode, TimeStep *tStep);

    virtual RemeshingCriteria *giveRemeshingCrit() { return NULL; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    const char *giveClassName() const { return "MeshQualityErrorEstimator"; }
    classType giveClassID() const { return MeshQualityErrorEstimatorClass; }
};
} // end namespace oofem
#endif // meshqualityerrorestimator_h






