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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef fecontactsurface_h
#define fecontactsurface_h

#include "contactsurface.h"
#include "contactprojection.h"
#include "contactelement.h"
#include "floatarray.h"
/*#include "chartype.h"
#include "domain.h"

#include "floatmatrix.h"
#include "feinterpol.h"
*/
//#include "valuemodetype.h"
//#include "unknowntype.h"
//#include "dofiditem.h"

/*
#include <cstdio>
#include <vector>
#include <memory>
*/
///@name Input fields for general element.
//@{
//@}

#define _IFT_FEContactSurface_contactElementSetNumber "ce_set"


namespace oofem {
class IntArray;
class ContactPoint;
/**
 * @brief Abstract representation of a finite element contact surface.
 *
 * The FEContactSurface class defines a common interface for contact surfaces
 * composed of finite elements. It encapsulates geometric and topological
 * information required for contact detection and enforcement, such as surface
 * discretization, local parameterization, and access to individual contact
 * points or segments.
 *
 * Concrete implementations represent specific surface discretizations
 * responsible for projecting a contact point onto a contact element contained
 * in the contact surface. This projection operation is a key component of
 * master–slave contact formulations and is used to determine closest-point
 * projections, gap functions, and local surface coordinates.
 *
 * FEContactSurface does not impose any particular contact enforcement method;
 * instead, it serves as an abstract geometric layer used by higher-level
 * contact algorithms and engineering models.
 */

class OOFEM_EXPORT FEContactSurface : public ContactSurface
{
protected:
    int contactElementSetNumber;
    /// Array containing dofmanager numbers.
    IntArray contactElementSet;
    /// Parametric-domain margin (in [-1,1]-normalized element coordinates)
    /// by which a projection may fall outside an element and still be
    /// accepted as "inside" it (comparable to other codes' sliding-elastic-
    /// interface search tolerance, typical default 0.01): a slave point
    /// whose ray/closest-point projection sits exactly on a shared edge
    /// between two facets will
    /// otherwise be rejected by one or both facets on alternating
    /// iterations, causing the owned facet to flip every Newton iteration.
    /// Default reproduces the previous unconditional near-zero tolerance.
    double domainTolerance = 1.e-10;

public:
    void setDomainTolerance(double tol) { domainTolerance = tol; }
    double giveDomainTolerance() const { return domainTolerance; }
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
    FEContactSurface(int n, Domain * aDomain);
    /// Virtual destructor.
    virtual ~FEContactSurface(){;}

    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;
    void postInitialize() override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
    /**
     * @brief Returns the i-th contact element associated with this surface.
     *
     * Provides access to a ContactElement by index in the internally stored list.
     * Indexing follows the internal ordering of the surface.
     *
     * @param i Index of the requested contact element.
     * @return Pointer to the contact element.
     */
    ContactElement* giveContactElement(int i);
    /**
     * @brief Returns the i-th contact element as referenced by the element set.
     *
     * Provides access to a ContactElement using the ordering of @c contactElementSet
     * (i.e., the element numbers stored from the configured set).
     *
     * @param i Index within the configured contact element set.
     * @return Pointer to the corresponding contact element.
 */
    ContactElement* giveContactElement_InSet(int i);
    /**
     * @brief Returns the number of contact elements on this surface.
     *
     * This corresponds to the size of the configured contact element set.
     *
     * @return Number of contact elements.
     */
    int giveNumberOfContactElements() const {return contactElementSet.giveSize();}
    /**
     * @brief Projects a contact point onto the smooth parameterization of a 3D contact element.
     *
     * Attempts to find an unconstrained stationary closest-point projection of @p cp onto the
     * smooth parameterization of @p e. A stationary projection on the closed parametric boundary
     * is admissible as a one-sided surface limit. An unconstrained projection outside that domain
     * is rejected: projecting it subsequently onto an edge or vertex would require separate
     * point-to-curve or point-to-corner kinematics and is not represented by the returned surface
     * normal and tangent basis.
     *
     * @param cp   Contact point to be projected (typically on the slave side).
     * @param e    Candidate master-side contact element on this surface.
     * @param tStep Current time step.
     *
     * @return Tuple with:
     *   - success flag (true if a valid, locally minimizing smooth-surface projection was found),
     *   - local coordinates (2D parametric coordinates on the element surface),
     *   - signed normal gap (convention defined by the contact formulation),
     *   - additional 3D vectors associated with the projection (e.g., normal and tangents).
     */
    virtual std::tuple <bool, FloatArrayF<2>, double, double, FloatArrayF<3>,FloatArrayF<3>,FloatArrayF<3>> findContactPointInElement_3d(ContactPoint *cp, ContactElement *e, TimeStep *tStep);
    /** Intersects a 3-D master facet with the ray through @p cp in @p direction. */
    virtual std::tuple <bool, FloatArrayF<2>, double, double, FloatArrayF<3>,FloatArrayF<3>,FloatArrayF<3>>
      findContactPointAlongDirectionInElement_3d(ContactPoint *cp, ContactElement *e,
                                                 const FloatArrayF<3> &direction,
                                                 TimeStep *tStep);
    /** Current unit normal of the finite-element surface carrying @p cp. */
    FloatArrayF<3> computeContactPointNormal_3d(ContactPoint *cp, TimeStep *tStep) const;
    /**
     * @brief Projects a contact point onto the smooth parameterization of a 2D contact element.
     *
     * An unconstrained stationary projection at an endpoint is admissible as a one-sided curve
     * limit. A projection whose unconstrained coordinate lies outside the segment is rejected;
     * endpoint contact outside the curve domain requires separate point-contact kinematics.
     *
     * @param cp   Contact point to be projected (typically on the slave side).
     * @param e    Candidate master-side contact element on this surface.
     * @param tStep Current time step.
     *
     * @return Tuple with:
     *   - success flag (true if a valid, locally minimizing smooth-curve projection was found),
     *   - local coordinate (1D parametric coordinate on the element),
     *   - signed normal gap (convention defined by the contact formulation),
     *   - additional 2D vectors associated with the projection (e.g., normal and tangent).
     */
    virtual std::tuple <bool, FloatArrayF<1>, double, double, FloatArrayF<2>,FloatArrayF<2>> findContactPointInElement_2d(ContactPoint *cp, ContactElement *e, TimeStep *tStep);
    /** Generalized point-to-edge projection for a linear 3-D contact facet. */
    ContactProjection findContactPointOnEdge_3d(ContactPoint *cp, ContactElement *e,
                                                int edgeIndex, TimeStep *tStep);
    /** Generalized point-to-vertex projection for a linear 3-D contact facet. */
    ContactProjection findContactPointOnVertex_3d(ContactPoint *cp, ContactElement *e,
                                                  int vertexIndex, TimeStep *tStep);
    /** Generalized point-to-corner projection for a linear 2-D contact segment. */
    ContactProjection findContactPointOnVertex_2d(ContactPoint *cp, ContactElement *e,
                                                  int vertexIndex, TimeStep *tStep);
    int giveNumberOfEdges(ContactElement *e) const;
    int giveNumberOfVertices(ContactElement *e) const;

private:
    /**
     * @brief Computes 3D local (parametric) coordinates of a projected contact point on a contact element.
     *
     * Internal helper performing the geometric projection of @p cp onto @p contactElement in 3D.
     * It returns whether the projection is admissible (e.g., inside the element parameter domain),
     * the local surface coordinates, the gap measure, and auxiliary geometric vectors used later
     * by the contact algorithm.
     *
     * @param cp             Contact point being projected.
     * @param contactElement Target contact element.
     * @param tStep          Current time step.
     *
     * @return Tuple with:
     *   - success flag,
     *   - local coordinates (2D parameter on the element surface),
     *   - signed normal gap,
     *   - additional 3D vectors associated with the projection (e.g., normal and tangent).
     */
    virtual std::tuple <bool, FloatArrayF<2>, double, double, FloatArrayF<3>,FloatArrayF<3>,FloatArrayF<3>> computeContactPointLocalCoordinates_3d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep);
    /**
     * @brief Computes 2D local (parametric) coordinate of a projected contact point on a contact element.
     *
     * Internal helper performing the geometric projection of @p cp onto @p contactElement in 2D.
     * It returns whether the projection is admissible (e.g., inside the element parameter domain),
     * the local coordinate, the gap measure, and auxiliary geometric vectors used later by the
     * contact algorithm.
     *
     * @param cp             Contact point being projected.
     * @param contactElement Target contact element.
     * @param tStep          Current time step.
     *
     * @return Tuple with:
     *   - success flag,
     *   - local coordinate (1D parameter on the element),
     *   - signed normal gap,
     *   - additional 2D vectors associated with the projection (e.g., normal and tangent).
     */
    virtual std::tuple <bool, FloatArrayF<1>, double, double, FloatArrayF<2>,FloatArrayF<2>>
      computeContactPointLocalCoordinates_2d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep);
    FloatArray computeIncidentPseudoNormal3d(ContactElement *representative,
                                             int firstNode, int secondNode,
                                             TimeStep *tStep) const;
    FloatArray computeIncidentPseudoNormal2d(ContactElement *representative,
                                             int vertexNode, TimeStep *tStep) const;

};


} // end namespace oofem
#endif //fecontactsurface_h
