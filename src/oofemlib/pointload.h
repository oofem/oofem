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
//   *** CLASS POINT LOAD ***
//   ************************

#ifndef pointload_h
#define pointload_h

#include "load.h"
#include "gausspnt.h"
#include "dictionr.h"

namespace oofem {

class TimeStep;

/**
 * Abstract base class representing a point load (force, momentum, ...) that acts
 * directly on or inside of some finite element.
 *
 * Methods for returning values of load components and returning the position of
 * load are provided.
 */
class PointLoad : public Load
{
    /*
     * This class implements a point load (force, moment,...) that acts
     * directly on or inside finite element.
     * A boundary load is usually attribute of one or more elements.
     * DESCRIPTION
     * The boundary load describes its geometry and values (it is assumed, that user will specify
     * all necessary dofs).
     *
     * REMARK
     * This class is not restricted to structural problems. For example, in ther-
     * mal analysis, a point load load could be a point heat source.
     */
public:

    // type to identify in which coordinate system
    // is load given
    /**
     * Load coordinate system type. Variable of this type can have following values BL_GlobalMode
     * (indicates that load given in global coordinate system) or BL_LocalMode
     * (entity dependent local coordinate system will be  used).
     */
    enum PL_CoordSystType {
        PL_GlobalMode, // global mode i.e. load is specifyied in global c.s.
        PL_LocalMode // local entity (edge or surface) coordinate system
    };

    /**
     * Type determining the type of formulation (entity local or global one).
     */
    enum PL_FormulationType {
        PL_EntityFormulation,
        PL_GlobalFormulation
    };

protected:
    /// Number of "DOFs" which represent load geometry
    int nDofs;
    /// Load type (its physical meaning)
    bcType lType;
    /// Load coordinate system
    PL_CoordSystType coordSystemType;
    /// aditional bc properties (coordinates, point of application)
    FloatArray coords;

public:
    /**
     * Constructor. Creates Boundary Load object with given number, belonging to given domain.
     * @param n load number
     * @param d domain to which new object will belongs.
     */
    PointLoad(int i, Domain *d) : Load(i, d) { nDofs = 0;
                                               coordSystemType = PL_GlobalMode; }                 // constructor

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
     * Returns coordinates of the reciver
     */
    void giveCoordinates(FloatArray &answer) { answer = coords; }
    /**
     * Return receiver's number of "DOFs". Should correspond to number of DOFs on loaded entity.
     */
    int                  giveNumberOfDofs() { return nDofs; }
    /**
     * Returns receiver's coordinate system
     */
    PL_CoordSystType     giveCoordSystMode() { return coordSystemType; }
    /**
     * Return formulation type
     */
    //virtual PL_FormulationType   giveFormulationType () {return BL_EntityFormulation;}
    /** Initializes receiver acording to object description stored in input record.
     *  Reads number of dofs into nDofs attribute (i.e. the number of dofs, which are on loaded entity),
     *  its loadType into loadType attribute and coordinate system type into csType attribute.
     */
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     *  @param str string to be filled by input record
     *  @param keyword print record keyword (default true)
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
    bcGeomType giveBCGeoType() const { return PointLoadBGT; }
    /** Returns classType id of receiver.
     *  @return BoundaryLoadClass value
     */
    classType            giveClassID() const { return PointLoadClass; }
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "PointLoad"; }

protected:
};

} // end namespace oofem
#endif // pointload_h
