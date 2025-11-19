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

#ifndef contactsurface_h
#define contactsurface_h

#include "femcmpnn.h"
#include "domain.h"

namespace oofem {
class IntArray;

/**
 * Abstract base class for all contact finite elements. Derived classes should be  base
 * classes for specific analysis type (for example base class for structural analysis,
 * thermal analysis or magnetostatics one). These derived classes then declare
 * analysis-specific part of interface and they provide default implementation
 * for these methods.
 * This abstract class declares (and possibly implements) general data and methods
 * common to all element types. General methods for obtaining characteristic vectors,
 * matrices and values are introduced and should be used instead of calling directly
 * specific member functions (these must be overloaded by derived analysis-specific
 * classes in order to invoke proper method according to type of component requested).
 */
class OOFEM_EXPORT ContactSurface : public FEMComponent
{
protected:
  
    /// Array containing dofmanager numbers.

public:
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
    ContactSurface(int n, Domain * aDomain);
    /// Virtual destructor.
    virtual ~ContactSurface(){;}

    /**@name Methods referring to code numbers */
    //@{
    /**
     * Returns the location array (array of code numbers) of receiv`er for given numbering scheme.
     * Results are cached at receiver for default scheme in locationArray attribute.
     */
    /*
    virtual void giveLocationArray(IntArray &locationArray, const UnknownNumberingScheme &s, IntArray *dofIds = NULL) const = 0;
    virtual void giveLocationArray(IntArray &locationArray, const IntArray &dofIDMask, const UnknownNumberingScheme &s, IntArray *dofIds = NULL) const = 0;
    */
    virtual void postInitialize() override {;}
    virtual void computeGaussPoints() { }

};


  
} // end namespace oofem
#endif //contactsurface_h
