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
class OOFEM_EXPORT ElementDofManager : public DofManager
{
private:
    Element *element;

public:
    /**
     * Constructor.
     * @param n Element dof manager number.
     * @param aDomain Domain which receiver belongs to.
     * @param elem Element to which receiver belongs.
     */
    ElementDofManager(int n, Domain * aDomain, Element * elem);
    /// Destructor.
    virtual ~ElementDofManager();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void printYourself();
    virtual const char *giveClassName() const { return "ElementDofManager"; }
    virtual const char *giveInputRecordName() const { return ""; } // Note: Can't be created in input files.

    virtual bool isDofTypeCompatible(dofType type) const { return ( type == DT_master || type == DT_simpleSlave ); }
};
} // end namespace oofem
#endif // elementdofman_h
