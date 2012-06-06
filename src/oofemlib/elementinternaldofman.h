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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef elementdofman_h
#define elementdofman_h

#include "dofmanager.h"

namespace oofem {

class Domain;
class Dof;
class NodalLoad;
class TimeStep;
class FloatArray;
class IntArray;
class Element;

/**
 * Class implementing internal element dof manager having some DOFs.
 * It possess degrees of freedom (see base class DofManager).
 * This class is usually attribute only of a single element, as its DOFs are internal element degrees of freedom.
 */
class ElementDofManager : public DofManager
{
private:
    Element* element;

public:
    /**
     * Constructor.
     * @param n Element dof manager number.
     * @param aDomain Domain which receiver belongs to.
     * @param elem Element to which receiver belongs.
     */
    ElementDofManager(int n, Domain *aDomain, Element* elem);
    /// Destructor.
    virtual ~ElementDofManager();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void printYourself();
    virtual const char *giveClassName() const { return "ElementDofManager"; }
    virtual classType giveClassID() const { return ElementDofManagerClass; }

    virtual bool isDofTypeCompatible(dofType type) const { return ( type == DT_master || type == DT_simpleSlave ); }
};
} // end namespace oofem
#endif // elementdofman_h






