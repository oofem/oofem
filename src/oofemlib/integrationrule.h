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

#ifndef integrationrule_h
#define integrationrule_h

#include "femcmpnn.h"
#include "materialmode.h"
#include "integrationdomain.h"
#include "contextioresulttype.h"
#include "contextmode.h"
#include "feinterpol.h"
#include "geometry.h"

namespace oofem {
class GaussPoint;
class Element;
class FEInterpolation;

/**
 * Abstract base class representing integration rule. The integration rule is
 * a collection of integration points used to  numerically integrate some formula.
 * The number of integration points and their coordinates and integration weights depends on
 * integration rule type (rule for integration in 1d, 2d, 3d) and required  accuracy.
 * General services for initialization are declared. Services for integration point retrieval are provided.
 *
 * In general, finite elements can have multiple integration rules, for different tasks or
 * when some components are integrated using reduced or selective integration.
 * Therefore, first and last index variables are introduced to characterize components
 * for which given integration rule applies.
 *
 * The integration rule is a rather passive object. It does not perform numerical
 * integration - it just provide way how to set up correct integration points
 * and weights.
 *
 * Because integration points contain related history parameters (using material status),
 * the unique copy of integration rule must exist on each element. The integration rule is
 * exclusively possessed by particular finite element.
 *
 *
 * Tasks:
 * - instanciating yourself
 * - returning number of integration points used
 * - returning requested integration point - method getIntegrationPoint
 * - returning interval of components (i.e.of local strain vector), where apply
 * - returning array of gauss points, according to specific
 *   integration rule (Gauss rule, Newton-Cotes rule, ...).
 *   integration points and corresponding weights are stored in Gauss point class.
 * - printing yourself
 * - updating yourself
 * - initializing for new time step
 * - saving & restoring context
 */
class IntegrationRule
{
protected:

    /// Number.
    int number;
    /// Element which integration rule is coupled to.
    Element *elem;

    /// Array containing integration points.
    GaussPoint **gaussPointArray;
    /// Number of integration point of receiver.
    int numberOfIntegrationPoints;
    /**
     * firstLocalStrainIndx and lastLocalStrainIndx indexes describe range of components (strains for example)
     * for which receiver integration points apply.
     */
    int firstLocalStrainIndx, lastLocalStrainIndx;

    /** 
     * Flag indicating that rule is dynamic, ie, its gauss points (their number, coordinates, weights) can change during
     * computation. Then some more data should be stored/restored from context file to reflect such dynamic feature.
     */
    bool isDynamic;

public:
    /**
     * Constructor.
     * @param n Number associated with receiver.
     * @param e Reference to element.
     * @param startIndx First component, for which rule applies.
     * @param endIndx Last component, for which rule applies.
     * @param dynamic Flag indicating that receiver can change.
     */
    IntegrationRule(int n, Element *e, int startIndx, int endIndx, bool dynamic);
    /**
     * Constructor.
     * @param n Number associated with receiver.
     * @param e Reference to element.
     */
    IntegrationRule(int n, Element *e);
    /// Destructor.
    virtual ~IntegrationRule();

    /**
     * Returns number of integration points of receiver.
     */
    int getNumberOfIntegrationPoints() const { return numberOfIntegrationPoints; }
    /**
     * Access particular integration point of receiver.
     * @param n Integration point number (should be in range 0,.., getNumberOfIntegrationPoints()-1).
     */
    GaussPoint *getIntegrationPoint(int n);
    /**
     * Returns starting component index, for which receiver applies.
     * @return First local strain index.
     */
    int getStartIndexOfLocalStrainWhereApply() { return firstLocalStrainIndx; }
    /**
     * Returns last component index, for which receiver applies.
     * @return Last local strain index.
     */
    int getEndIndexOfLocalStrainWhereApply() { return lastLocalStrainIndx; }
    /**
     * Initializes the receiver. Receiver integration points are created according to given parameters.
     * @param mode Describes integration domain.
     * @param nPoints Required number of integration points of receiver.
     * @param matMode Material mode of receiver's integration points.
     * @return Number of points.
     */
    int setUpIntegrationPoints(integrationDomain mode, int nPoints, MaterialMode matMode);
    /**
     * Initializes the receiver. Receiver integration points are created according to given parameters.
     * @param mode Describes integration domain.
     * @param nPoints Required number of integration points of receiver.
     * @param matMode Material mode of receiver's integration points.
     * @param coords
     * @return Number of points.
     */
    int setUpEmbeddedIntegrationPoints(integrationDomain mode, int nPoints, MaterialMode matMode,
                                       const FloatArray **coords);

    /**
     * Prints receiver's output to given stream.
     * Invokes printOutputAt service on all receiver's integration points.
     */
    virtual void printOutputAt(FILE *file, TimeStep *stepN);
    /**
     * Updates receiver state.
     * Calls updateYourself service of all receiver's integration points.
     */
    void updateYourself(TimeStep *tStep);
    /**
     * Initializes receiver.
     * Calls initForNewStep service of all receiver's integration points.
     */
    void initForNewStep();

    /** Returns reference to element containing receiver */
    Element *giveElement() { return elem; }
    /** Returns reference to interpolation associated with GaussPoint */
    FEInterpolation *giveInterpolation(GaussPoint *gp);
    /** Returns receiver number */
    int giveNumber() { return this->number; }
    /**
     * Abstract service.
     * Returns required number of integration points to exactly integrate
     * polynomial of order approxOrder on given domain.
     * When approxOrder is too large and is not supported by implementation
     * method returns -1. Must be overloaded by derived classes.
     */
    virtual int getRequiredNumberOfIntegrationPoints(integrationDomain dType, int approxOrder) { return 0; }

    /**
     * Saves receiver's context to stream.
     * Calls saveContext service for all receiver's integration points.
     * Note: does not call the FEMComponent::saveContext service, in order not
     * to write class id info for each integration rule.
     * @param stream Output stream.
     * @param mode Determines amount of info required in stream (state, definition, ...).
     * @param obj Special parameter.
     * @exception ContextIOERR If error encountered.
     */
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj);
    /**
     * Restores receiver's context to stream.
     * Calls restoreContext service for all receiver's integration points.
     * Note: does not call the FEMComponent::restoreContext service, in order not
     * to write class id info for each integration rule.
     * @param stream Input stream.
     * @param mode Determines amount of info available in stream (state, definition, ...).
     * @param obj Should be a pointer to invoking element, ie., to which the receiver will belong to.
     * @exception ContextIOERR If error encountered.
     */
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj);
    /**
     * Clears the receiver, ie deallocates all integration points
     */
    void clear();

    /// Returns receiver sub patch indices (if apply).
    virtual const IntArray *giveKnotSpan() { return NULL; }

    /// Returns classType id of receiver.
    virtual classType giveClassID() const { return IntegrationRuleClass; }
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "IntegrationRule"; }
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }

protected:
    /**
     * Sets up receiver's  integration points on unit line integration domain.
     * Default implementation does not sets up any integration points and returns 0.
     * Must be overloaded by derived classes.
     * @returns Number of integration points.
     */
    virtual int SetUpPointsOnLine(int, MaterialMode mode, GaussPoint ***gp) { return 0; }
    /**
     * Sets up receiver's  integration points on triangular (area coords) integration domain.
     * Default implementation does not sets up any integration points and returns 0.
     * Must be overloaded by derived classes.
     * @returns Number of integration points.
     */
    virtual int SetUpPointsOnTriagle(int, MaterialMode mode, GaussPoint ***gp) { return 0; }
    /**
     * Sets up receiver's  integration points on unit square integration domain.
     * Default implementation does not sets up any integration points and returns 0.
     * Must be overloaded by derived classes.
     * @returns Number of integration points.
     */
    virtual int SetUpPointsOnSquare(int, MaterialMode mode, GaussPoint ***gp) { return 0; }
    /**
     * Sets up receiver's  integration points on unit cube integration domain.
     * Default implementation does not sets up any integration points and returns 0.
     * Must be overloaded by derived classes.
     * @returns Number of integration points.
     */
    virtual int SetUpPointsOnCube(int, MaterialMode mode, GaussPoint ***gp) { return 0; }
    /**
     * Sets up receiver's  integration points on tetrahedra (volume coords) integration domain.
     * Default implementation does not sets up any integration points and returns 0.
     * Must be overloaded by derived classes.
     * @returns Number of integration points.
     */
    virtual int SetUpPointsOnTetrahedra(int, MaterialMode mode, GaussPoint ***gp) { return 0; }
    /**
     * Sets up integration points on 2D embedded line inside 2D volume (the list of local coordinates
     * should be provided).
     */
    virtual int SetUpPointsOn2DEmbeddedLine(int nPoints, MaterialMode mode, GaussPoint ***,
                                            const FloatArray **coords) { return 0; }
};
} // end namespace oofem
#endif // integrationrule_h
