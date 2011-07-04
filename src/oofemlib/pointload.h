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
 *
 * @note
 * This class is not restricted to structural problems. For example, in thermal
 * analysis, a point load load could be a point heat source.
 */
class PointLoad : public Load
{
public:
    /**
     * Load coordinate system type. Variable of this type can have following values BL_GlobalMode
     * (indicates that load given in global coordinate system) or BL_LocalMode
     * (entity dependent local coordinate system will be  used).
     */
    enum PL_CoordSystType {
        PL_GlobalMode, ///< Global mode i.e. load is specified in global c.s.
        PL_LocalMode, ///< Local entity (edge or surface) coordinate system.
    };

    /**
     * Type determining the type of formulation (entity local or global one).
     */
    enum PL_FormulationType {
        PL_EntityFormulation,
        PL_GlobalFormulation,
    };

protected:
    /// Number of "DOFs" which represent load geometry.
    int nDofs;
    /// Load type (its physical meaning).
    bcType lType;
    /// Load coordinate system.
    PL_CoordSystType coordSystemType;
    /// Additional properties (coordinates, point of application).
    FloatArray coords;

public:
    /**
     * Constructor. Creates Boundary Load object with given number, belonging to given domain.
     * @param n Load number.
     * @param d Domain to which new object will belongs.
     */
    PointLoad(int n, Domain *d) : Load(n, d) {
        nDofs = 0;
        coordSystemType = PL_GlobalMode;
    }

    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, FloatArray &coords, ValueModeType mode);
    /**
     * Gives coordinates of the receiver
     */
    void giveCoordinates(FloatArray &answer) { answer = coords; }
    /**
     * Return receiver's number of "DOFs". Should correspond to number of DOFs on loaded entity.
     */
    int giveNumberOfDofs() { return nDofs; }
    /**
     * Returns receiver's coordinate system
     */
    PL_CoordSystType giveCoordSystMode() { return coordSystemType; }
    /*
     * Return formulation type
     */
    //virtual PL_FormulationType giveFormulationType () {return BL_EntityFormulation;}

    IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    bcType giveType() const { return lType; }
    bcGeomType giveBCGeoType() const { return PointLoadBGT; }
    classType giveClassID() const { return PointLoadClass; }
    const char *giveClassName() const { return "PointLoad"; }

protected:
};
} // end namespace oofem
#endif // pointload_h
