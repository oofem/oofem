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
    ContactElement* giveContactElement(int i);
    ContactElement* giveContactElement_InSet(int i);
    int giveNumberOfContactElements(){return contactElementSet.giveSize();}
    virtual std::tuple <bool, FloatArrayF<2>, double,FloatArrayF<3>,FloatArrayF<3>,FloatArrayF<3>> findContactPointInElement_3d(ContactPoint *cp, ContactElement *e, TimeStep *tStep);
    virtual std::tuple <bool, FloatArrayF<1>, double,FloatArrayF<2>,FloatArrayF<2>> findContactPointInElement_2d(ContactPoint *cp, ContactElement *e, TimeStep *tStep);

private:
    virtual std::tuple <bool, FloatArrayF<2>, double,  FloatArrayF<3>,FloatArrayF<3>,FloatArrayF<3>> computeContactPointLocalCoordinates_3d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep);
    virtual std::tuple <bool,  FloatArrayF<1>, double, FloatArrayF<2>,FloatArrayF<2>> 
      computeContactPointLocalCoordinates_2d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep);

};


} // end namespace oofem
#endif //contactelement_h
