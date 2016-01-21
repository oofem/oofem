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

#ifndef integrationrule_h
#define integrationrule_h

#include "oofemcfg.h"
#include "materialmode.h"
#include "integrationdomain.h"
#include "contextioresulttype.h"
#include "contextmode.h"
#include "inputrecord.h"

#include <cstdio>

namespace oofem {
class TimeStep;
class GaussPoint;
class Element;
class DataStream;

///@todo Breaks modularity, reconsider this;
enum IntegrationRuleType {
    IRT_None = 0,
    IRT_Gauss = 1,
    IRT_Lobatto = 2
};

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
class OOFEM_EXPORT IntegrationRule
{
protected:
    /// Number.
    int number;
    /// Element which integration rule is coupled to.
    Element *elem;
    /// Integration domain
    integrationDomain intdomain;

    /// Array containing integration points.
    std::vector< GaussPoint * > gaussPoints;
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
    std::vector< GaussPoint *> :: iterator begin() { return gaussPoints.begin(); }
    std::vector< GaussPoint *> :: iterator end() { return gaussPoints.end(); }

    /**
     * Constructor.
     * @param n Number associated with receiver.
     * @param e Reference to element.
     * @param startIndx First component, for which rule applies.
     * @param endIndx Last component, for which rule applies.
     * @param dynamic Flag indicating that receiver can change.
     */
    IntegrationRule(int n, Element * e, int startIndx, int endIndx, bool dynamic);
    /**
     * Constructor.
     * @param n Number associated with receiver.
     * @param e Reference to element.
     */
    IntegrationRule(int n, Element * e);
    /// Destructor.
    virtual ~IntegrationRule();

    /**
     * Returns number of integration points of receiver.
     */
    int giveNumberOfIntegrationPoints() const { return (int)gaussPoints.size(); }
    /**
     * Access particular integration point of receiver.
     * @param n Integration point number (should be in range 0,.., giveNumberOfIntegrationPoints()-1).
     */
    GaussPoint *getIntegrationPoint(int n);
    /**
     * Scans through the integration points and finds the one closest to the given (local) coordinate.
     */
    GaussPoint *findIntegrationPointClosestTo(const FloatArray &lcoord);
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
     * @param intdomain Describes integration domain.
     * @param nPoints Required number of integration points of receiver.
     * @param matMode Material mode of receiver's integration points.
     * @return Number of points.
     */
    int setUpIntegrationPoints(integrationDomain intdomain, int nPoints, MaterialMode matMode);
    /**
     * Initializes the receiver. Receiver integration points are created according to given parameters.
     * @param intdomain Describes integration domain.
     * @param nPoints Required number of integration points of receiver.
     * @param matMode Material mode of receiver's integration points.
     * @param coords
     * @return Number of points.
     */
    int setUpEmbeddedIntegrationPoints(integrationDomain intdomain, int nPoints, MaterialMode matMode,
                                       const std :: vector< FloatArray > &coords);

    /**
     * Prints receiver's output to given stream.
     * Invokes printOutputAt service on all receiver's integration points.
     */
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    /**
     * Updates receiver state.
     * Calls updateYourself service of all receiver's integration points.
     */
    void updateYourself(TimeStep *tStep);

    /** Returns reference to element containing receiver */
    Element *giveElement() { return elem; }
    /** Returns receiver number */
    int giveNumber() { return this->number; }
    /** Returns the domain for the receiver */
    integrationDomain giveIntegrationDomain() const { return this->intdomain; }
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
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj);
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
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj);
    /**
     * Clears the receiver, ie deallocates all integration points
     */
    void clear();

    /// Returns receiver sub patch indices (if apply).
    virtual const IntArray *giveKnotSpan() { return NULL; }

    virtual const char *giveClassName() const { return "IntegrationRule"; }
    /// Error printing helper.
    std :: string errorInfo(const char *func) const { return std :: string(giveClassName()) + func; }
    virtual IntegrationRuleType giveIntegrationRuleType() const { return IRT_None; }
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }

    /**
     * Trivial implementation, only creates a single point.
     */
    int SetUpPoint(MaterialMode mode);
    /**
     * Sets up receiver's integration points on unit line integration domain.
     * Default implementation does not sets up any integration points and returns 0.
     * Must be overloaded by derived classes.
     * @return Number of integration points.
     */
    virtual int SetUpPointsOnLine(int, MaterialMode mode) { return 0; }
    /**
     * Sets up receiver's integration points on triangular (area coords) integration domain.
     * Default implementation does not sets up any integration points and returns 0.
     * Must be overloaded by derived classes.
     * @return Number of integration points.
     */
    virtual int SetUpPointsOnTriangle(int, MaterialMode mode) { return 0; }
    /**
     * Sets up receiver's integration points on unit square integration domain.
     * Default implementation does not sets up any integration points and returns 0.
     * Must be overloaded by derived classes.
     * @return Number of integration points.
     */
    virtual int SetUpPointsOnSquare(int, MaterialMode mode) { return 0; }
    /**
     * Sets up receiver's integration points on unit cube integration domain.
     * Default implementation does not sets up any integration points and returns 0.
     * Must be overloaded by derived classes.
     * @return Number of integration points.
     */
    virtual int SetUpPointsOnCube(int, MaterialMode mode) { return 0; }
    /**
     * Sets up receiver's integration points on unit cube integration domain divided into layers in the zeta-direction.
     * Default implementation does not sets up any integration points and returns 0.
     * Must be overloaded by derived classes.
     * @param nPoints1 Number of integration points in the "xi"-direction.
     * @param nPoints2 Number of integration points in the "eta"-direction.
     * @param nPointsDepth Number of integration points in the "zeta"-direction
     * @return Number of integration points.
     */
    virtual int SetUpPointsOnCubeLayers(int nPoints1, int nPoints2, int nPointsDepth, MaterialMode mode, const FloatArray &layerThickness) { return 0; }
    /**
     * Sets up receiver's integration points on tetrahedra (volume coords) integration domain.
     * Default implementation does not sets up any integration points and returns 0.
     * Must be overloaded by derived classes.
     * @return Number of integration points.
     */
    virtual int SetUpPointsOnTetrahedra(int, MaterialMode mode) { return 0; }
    /**
     * Sets up integration points on 2D embedded line inside 2D volume (the list of local coordinates
     * should be provided).
     * @param nPoints Number of points along line.
     */
    virtual int SetUpPointsOn2DEmbeddedLine(int nPoints, MaterialMode mode,
                                            const FloatArray &coord0, const FloatArray &coord1) { return 0; }

    /**
     * Sets up receiver's integration points on a wedge integration domain.
     * Default implementation does not sets up any integration points and returns 0.
     * Must be overloaded by derived classes.
     * @param nPointsTri Number of points over the triangle cross-section.
     * @param nPointsDepth Number of points over the depth.
     * @return Number of integration points.
     */
    virtual int SetUpPointsOnWedge(int nPointsTri, int nPointsDepth, MaterialMode mode) { return 0; }
    /**
     * Sets up receiver's integration points on a wedge integration domain divided into layers in the zeta-direction.
     * Default implementation does not sets up any integration points and returns 0.
     * Must be overloaded by derived classes.
     * @param nPointsTri Number of points over the triangle cross-section.
     * @param nPointsDepth Number of points over the depth.
     * @return Number of integration points.
     */
    virtual int SetUpPointsOnWedgeLayers(int nPointsTri, int nPointsDepth, MaterialMode mode, const FloatArray &layerThickness) { return 0; }
};
} // end namespace oofem
#endif // integrationrule_h
