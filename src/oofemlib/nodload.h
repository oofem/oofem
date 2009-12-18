/* $Header: /home/cvs/bp/oofem/oofemlib/src/nodload.h,v 1.6 2003/04/06 14:08:25 bp Exp $ */
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
//   *** CLASS NODAL LOAD ***
//   ************************

#ifndef nodload_h
#define nodload_h

#include "load.h"
#include "bcgeomtype.h"
#include "valuemodetype.h"

namespace oofem {

class TimeStep;

/**
 * Class implementing a concentrated load (force, moment,...) that acts
 * directly on a dof manager (node or element side, if it has associated DOFs).
 * This load could not be applied on an element.
 * A nodal load is usually attribute of one or more nodes or element sides.
 *
 * The component array, which size should be same as number of DOFs in particular
 * node/side is read from input.
 */
class NodalLoad : public Load
{
    /*
     * This class implements a concentrated load (force, moment,...) that acts
     * directly on a node (not on an element). A nodal load is usually attribute
     * of one or more nodes.
     * DESCRIPTION
     * The attribute 'componentArray' contains, for example for the case of a
     * plane beam structure, 2 forces and 1 moment on a right place. (6 dof per
     * node is assumed)
     * REMARK
     * class Node is not restricted to structural problems. For example, in ther-
     * mal analysis, a nodal load would be a concentrated heat source.
     */
public:
    /**
     * Type determining the type of formulation (entity local or global one).
     */
    enum BL_CoordSystType {
        BL_GlobalMode, // global mode i.e. load is specifyied in global c.s.
        BL_LocalMode // local entity (edge or surface) coordinate system
    };

protected:
    /** Load coordinate system
     *  it is actually used only when local coordinate system in node is defined and load is specified in global
     *  coordinate system */
    BL_CoordSystType coordSystemType;



public:

    /**
     * Constructor. Creates nodal load object with given number, belonging to given domain.
     * @param n load  number
     * @param d domain to which new object will belongs.
     */
    NodalLoad(int i, Domain *d) : Load(i, d) { }      // constructor
    /**
     * Returns receiver bc type. It distinguish particular boundary conditions according to
     * their "physical" meaning (like StructuralTemperatureLoadLT, StructuralLoadLT).
     * @return returns DofManagerLoadLT value.
     */
    bcGeomType     giveBCGeoType() const { return NodalLoadBGT; }
    /// Returns input record name of the receiver.
    const char *giveInputRecordName() const { return "NodalLoad"; }
    void         computeValueAt(FloatArray &answer, TimeStep *atTime, FloatArray &coords, ValueModeType mode)
    { computeComponentArrayAt(answer, atTime, mode); }
    /**
     * Returns receiver's coordinate system
     */
    BL_CoordSystType     giveCoordSystMode() { return coordSystemType; }
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
};

} // end namespace oofem
#endif // nodload_h
