/* $Header: /home/cvs/bp/oofem/oofemlib/src/boundaryload.h,v 1.10 2003/04/06 14:08:23 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   ************************
//   *** CLASS BOUNDARY LOAD ***
//   ************************

#ifndef boundaryload_h
#define boundaryload_h

#include "load.h"
#include "gausspnt.h"
#include "dictionr.h"

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
 * then its values are evalueated at points, which is characterized by global coordinates.
 * If load is specified in entity space, then point is characterized by entity isoparametric coordinates.
 *
 * Methods for evaluation of load component values in any point on element boundary (on side, face, ...)
 * are provided (this point is determined using global or local (global or entity definition) coordinates).
 * The other possibility (faster, but less general) is to specify the the point using isoparametric
 * coordinates - but this is not supported.
 * It is generally assumed, that derived classes will approximate somehow their values
 * based on user specified data (on side nodes for example) and load approximation type.
 * The similar scheme borrowed from FE appriximation is used, values computed at required point
 * are computed as a product of approximation matrix (matrix of approximation functions) with
 * "vertex" values, which has to be specified on input by user.
 * Elements can request the order of load approximation (for setting up the appropriate
 * integration rule order) and the array of component values (for each dof) at specific integration point
 * on the boundary.
 *
 * For some elements it may be better to obtain "vertex values" of boundary load to
 * compute load vector directly using exact formulae. Elements then can ask for
 * values at nodal points and obtain cooresponding  "vertex values". Meaning of these values is
 * class dependent, see derived classes documentation for details.
 *
 * Elements must take care, on which boundary the load acts on (side number, ...).
 * Boundary load class also introduces load related cooordinate system indicator.
 * Load can be generally specified in global coordinate system or in entity dependent local coordinate
 * system. The entity dependent coordinate system is defined by particular element.
 *
 * Note, this class is not restricted to structural problems. For example, in ther-
 * mal analysis, a boundary load load could be a  heat source.
 *
 * To sumarize, the services provided include
 * <UL>
 * <LI>
 * Computing component array evaluated at specific point on boundary.</LI>
 * <LI>
 * Returning component array of "vertex values". Meaning of these values is
 * class dependent, see derived classes documentation for details. "vertexes" can generally differ from
 * element nodes.</LI>
 * <LI>
 * Returning load appriximation order. Usefull when numerical integrations of load vector over element
 * boundaries are used.</LI>
 * <LI>
 * Returning type of coordinate system, in which load applies. (global c.s., or entity related c.s.).</LI>
 * </UL>
 * Returning number of load approximation DOFs, which represent its geometry.
 * (number of DOFs is also size of load component array attribute and should correspond to number of DOFs on loaded
 * entity).
 */
class BoundaryLoad : public Load
{
    /*
     * This class implements a boundary load (force, moment,...) that acts
     * directly on a boundary of some finite element (on side, face, ..).
     * A boundary load is usually attribute of one or more elements.
     * DESCRIPTION
     * The boundary load describes its geometry and values (it is assumed, that user will specify
     * all necessary dofs) on  boundary using isoparametric approximation.
     * Elements can request the order of approximation (for setting up the appropriate
     * integration rule order) and the array of values (for each dof) at specific integration point
     * on the boundary.
     *
     * Elements must take care, on which boundary the load acts on (side number, ...).
     *
     * For some elements it may be better to obtain "vertex values" of boundary load to
     * compute load vector directly using exact formulae.
     *
     *
     * REMARK
     * This class is not restricted to structural problems. For example, in ther-
     * mal analysis, a boundary load load would be a  heat source.
     */
public:

    // type to identify in which coordinate system
    // is load given
    /**
     * Load coordinate system type. Variable of this type can have following values BL_GlobalMode
     * (indicates that load given in global coordinate system) or BL_LocalMode
     * (entity dependent local coordinate system will be  used).
     */
    enum BL_CoordSystType {
        BL_GlobalMode, // global mode i.e. load is specifyied in global c.s.
        BL_LocalMode // local entity (edge or surface) coordinate system
    };

    /**
     * Type determining the type of formulation (entity local or global one).
     */
    enum BL_FormulationType {
        BL_EntityFormulation,
        BL_GlobalFormulation
    };

protected:
    /// Number of "DOFs" which represent load geometry
    int nDofs;
    /// Load type (its physical meaning)
    bcType lType;
    /// Load coordinate system
    BL_CoordSystType coordSystemType;
    /// aditional bc properties
    Dictionary propertyDictionary;

public:
    /**
     * Constructor. Creates Boundary Load object with given number, belonging to given domain.
     * @param n load number
     * @param d domain to which new object will belongs.
     */
    BoundaryLoad(int i, Domain *d) : Load(i, d) { nDofs = 0;
                                                  coordSystemType = BL_GlobalMode; }                 // constructor

    /**
     * Computes components values of load at given point - global coordinates (coordinates given).
     * Default implementation computes product of aproximation matrix (computeNArray service) and
     * with "vertex" value array attribute and the result is then multiplied by
     * corresponding load time function value respecting load response mode.
     * @param answer component values at given point and time
     * @param stepN time step representing time
     * @param coords global (or local) problem coordinates, which are used to
     * evaluate components values.
     * @param mode determines response mode.
     */
    virtual void         computeValueAt(FloatArray &answer, TimeStep *, FloatArray &coords, ValueModeType mode);
    /**
     * Returns approximation order of load geometry.
     */
    virtual int          giveApproxOrder() = 0;
    /**
     * Return receiver's number of "DOFs". Should correspond to number of DOFs on loaded entity.
     */
    int                  giveNumberOfDofs() { return nDofs; }
    /**
     * Returns receiver's coordinate system
     */
    BL_CoordSystType     giveCoordSystMode() { return coordSystemType; }
    /**
     * Return formulation type
     */
    virtual BL_FormulationType   giveFormulationType() { return BL_EntityFormulation; }
    /** Initializes receiver acording to object description stored in input record.
     * Reads number of dofs into nDofs attribute (i.e. the number of dofs, which are on loaded entity),
     * its loadType into loadType attribute and coordinate system type into csType attribute.
     */
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /**
     * Returns receiver load type. It distinguish particular boundary conditions according to
     * their "physical" meaning (like StructuralTemperatureLoadLT, StructuralLoadLT).
     * Derived classes should always overload, default implementation returns value
     * specified on input by user.
     * See cltypes.h file for details.
     */
    bcType               giveType() const { return lType; }
    /**
     * Returns the value of a property 'aProperty'. Property must be identified
     * by unique int id.
     * @aProperty id of peroperty requested
     * @return property value
     */
    virtual double   giveProperty(int);

    /** Returns classType id of receiver.
     * @return BoundaryLoadClass value
     */
    classType            giveClassID() const { return BoundaryLoadClass; }
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "BoundaryLoad"; }

protected:
    /**
     * Abstract function, for computing approximation matrix of receiver at given  point.
     * The product of approximation matrix with "vertex" values array atrribute will produce
     * load components in given  point.
     * @param answer approximation matrix
     * @param coords global integration point coordinates
     */
    virtual void computeNArray(FloatArray &answer, FloatArray &coords) const  = 0;
    /**
     * Returns array of load "vertex" values evaluated at given time.
     * @param answer load "vertex" values
     * @param tStep time step
     * @param mode determines response mode
     */
    void                 computeComponentArrayAt(FloatArray &answer, TimeStep *, ValueModeType mode);
};

} // end namespace oofem
#endif // boundaryload_h






