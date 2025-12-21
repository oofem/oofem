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

public:
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
    FEContactSurface(int n, Domain * aDomain);
    /// Virtual destructor.
    virtual ~FEContactSurface(){;}

    void initializeFrom(InputRecord &ir) override;
    void postInitialize() override;
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
    int giveNumberOfContactElements(){return contactElementSet.giveSize();}
    /**
     * @brief Projects a contact point onto a 3D contact element and evaluates contact kinematics.
     *
     * Attempts to find the projection of the contact point @p cp onto the 3D contact element @p e
     * for the given time step. On success, returns the local surface coordinates of the projection,
     * the normal gap measure, and associated geometric vectors required by the contact formulation
     * (e.g., normal/tangent directions depending on the element type).
     *
     * @param cp   Contact point to be projected (typically on the slave side).
     * @param e    Candidate master-side contact element on this surface.
     * @param tStep Current time step.
     *
     * @return Tuple with:
     *   - success flag (true if a valid projection was found),
     *   - local coordinates (2D parametric coordinates on the element surface),
     *   - signed normal gap (convention defined by the contact formulation),
     *   - additional 3D vectors associated with the projection (e.g., normal and tangents).
     */
    virtual std::tuple <bool, FloatArrayF<2>, double,FloatArrayF<3>,FloatArrayF<3>,FloatArrayF<3>> findContactPointInElement_3d(ContactPoint *cp, ContactElement *e, TimeStep *tStep);
    /**
     * @brief Projects a contact point onto a 2D contact element and evaluates contact kinematics.
     *
     * Attempts to find the projection of the contact point @p cp onto the 2D contact element @p e
     * for the given time step. On success, returns the local coordinate of the projection (1D
     * parameter), the normal gap measure, and associated geometric vectors required by the contact
     * formulation (e.g., normal/tangent directions depending on the element type).
     *
     * @param cp   Contact point to be projected (typically on the slave side).
     * @param e    Candidate master-side contact element on this surface.
     * @param tStep Current time step.
     *
     * @return Tuple with:
     *   - success flag (true if a valid projection was found),
     *   - local coordinate (1D parametric coordinate on the element),
     *   - signed normal gap (convention defined by the contact formulation),
     *   - additional 2D vectors associated with the projection (e.g., normal and tangent).
     */
    virtual std::tuple <bool, FloatArrayF<1>, double,FloatArrayF<2>,FloatArrayF<2>> findContactPointInElement_2d(ContactPoint *cp, ContactElement *e, TimeStep *tStep);

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
    virtual std::tuple <bool, FloatArrayF<2>, double,  FloatArrayF<3>,FloatArrayF<3>,FloatArrayF<3>> computeContactPointLocalCoordinates_3d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep);
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
    virtual std::tuple <bool,  FloatArrayF<1>, double, FloatArrayF<2>,FloatArrayF<2>> 
      computeContactPointLocalCoordinates_2d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep);

};


} // end namespace oofem
#endif //fecontactsurface_h
