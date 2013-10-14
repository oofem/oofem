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

#ifndef boundaryload_h
#define boundaryload_h

#include "load.h"
#include "bctype.h"
#include "valuemodetype.h"
#include "dictionary.h"

///@name Input fields for BoundaryLoad
//@{
#define _IFT_BoundaryLoad_ndofs "ndofs"
#define _IFT_BoundaryLoad_loadtype "loadtype"
#define _IFT_BoundaryLoad_cstype "cstype"
#define _IFT_BoundaryLoad_properties "properties"
//@}

namespace oofem {
class TimeStep;

/**
 * Abstract base class representing a boundary load (force, momentum, ...) that acts
 * directly on a boundary of some finite element (on element side, face,  ...).
 * Boundary load is usually attribute of one or more finite elements.
 *
 * This base class only declares the common services common to all derived classes.
 * Derived classes must implement abstract services and possibly may customize existing.
 * Boundary load is represented by its geometry (determined by its type - linear, quadratic load)
 * and values (it is assumed, that user will supply all necessary values for each dof).
 *
 * The load can generally be specified in global space or can be related
 * to local entity space (related to edge, surface). If load is specified in global space
 * then its values are evaluated at points, which is characterized by global coordinates.
 * If load is specified in entity space, then point is characterized by entity isoparametric coordinates.
 *
 * Methods for evaluation of load component values in any point on element boundary (on side, face, ...)
 * are provided (this point is determined using global or local (global or entity definition) coordinates).
 * The other possibility (faster, but less general) is to specify the the point using isoparametric
 * coordinates - but this is not supported.
 * It is generally assumed, that derived classes will approximate somehow their values
 * based on user specified data (on side nodes for example) and load approximation type.
 * The similar scheme borrowed from FE approximation is used, values computed at required point
 * are computed as a product of approximation matrix (matrix of approximation functions) with
 * "vertex" values, which has to be specified on input by user.
 * Elements can request the order of load approximation (for setting up the appropriate
 * integration rule order) and the array of component values (for each dof) at specific integration point
 * on the boundary.
 *
 * For some elements it may be better to obtain "vertex values" of boundary load to
 * compute load vector directly using exact formulae. Elements then can ask for
 * values at nodal points and obtain corresponding  "vertex values". Meaning of these values is
 * class dependent, see derived classes documentation for details.
 *
 * Elements must take care, on which boundary the load acts on (side number, ...).
 * Boundary load class also introduces load related coordinate system indicator.
 * Load can be generally specified in global coordinate system or in entity dependent local coordinate
 * system. The entity dependent coordinate system is defined by particular element.
 *
 * Note, this class is not restricted to structural problems. For example, in thermal
 * analysis, a boundary load load could be a  heat source.
 *
 * To summarize, the services provided include
 * - Computing component array evaluated at specific point on boundary.
 * - Returning component array of "vertex values". Meaning of these values is
 *   class dependent, see derived classes documentation for details. "vertexes" can generally differ from
 *   element nodes.
 * - Returning load approximation order. Useful when numerical integrations of load vector over element
 *   boundaries are used.
 * - Returning type of coordinate system, in which load applies. (global c.s., or entity related c.s.).
 * - Returning number of load approximation DOFs, which represent its geometry.
 *   (number of DOFs is also size of load component array attribute and should correspond to number of DOFs on loaded
 *   entity).
 */
class OOFEM_EXPORT BoundaryLoad : public Load
{
protected:
    /// Number of "DOFs" which represent load geometry.
    int nDofs;
    /// Load type (its physical meaning).
    bcType lType;
    /// Load coordinate system.
    CoordSystType coordSystemType;
    /// Additional b.c properties.
    Dictionary propertyDictionary;

public:
    /**
     * Constructor. Creates a boundary load object with given number, belonging to given domain.
     * @param i Load number.
     * @param d Domain to which new object will belongs.
     */
    BoundaryLoad(int i, Domain *d) : Load(i, d), nDofs(0), coordSystemType(CST_Global) { }

    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, FloatArray &coords, ValueModeType mode);
    /**
     * @return Approximation order of load geometry.
     */
    virtual int giveApproxOrder() = 0;
    /**
     * @return Receiver's number of "DOFs". Should correspond to number of DOFs on loaded entity.
     */
    int giveNumberOfDofs() { return nDofs; }

    virtual CoordSystType giveCoordSystMode() { return coordSystemType; }
    /**
     * Initializes receiver according to object description stored in input record.
     * Reads number of dofs into nDofs attribute (i.e. the number of dofs, which are on loaded entity),
     * its loadType into loadType attribute and coordinate system type into csType attribute.
     */
    IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    /**
     * Returns receiver load type. It distinguish particular boundary conditions according to
     * their "physical" meaning (like StructuralTemperatureLoadLT, StructuralLoadLT).
     * Derived classes should always overload, default implementation returns value
     * specified on input by user.
     * See cltypes.h file for details.
     */
    virtual bcType giveType() const { return lType; }
    virtual double giveProperty(int aProperty);

protected:
    /**
     * Abstract function, for computing approximation matrix of receiver at given  point.
     * The product of approximation matrix with "vertex" values array attribute will produce
     * load components in given  point.
     * @param answer Approximation vector.
     * @param coords Global integration point coordinates.
     */
    virtual void computeNArray(FloatArray &answer, FloatArray &coords) const  = 0;
    /**
     * Returns array of load "vertex" values evaluated at given time.
     * @param answer Load "vertex" values.
     * @param tStep Time step.
     * @param mode Determines response mode.
     */
    virtual void computeComponentArrayAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
};
} // end namespace oofem
#endif // boundaryload_h






