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


#ifndef latticetransportelement_h
#define latticetransportelement_h

#include "transportelement.h"

/**
 * This class implements the base of a special transport lattice element following
 * the concepts orginally developed by John Bolander. In this lattice
 * framework, elements are pipes used for either heat transfer or mass transport.
 * In this base class common interfaces of derived elements are defined.
 */

namespace oofem {
class LatticeTransportElement : public TransportElement
{
public:
    LatticeTransportElement(int, Domain *, ElementMode);        // constructor

    ~LatticeTransportElement();                                 // destructor

    IRResultType initializeFrom(InputRecord *ir);

    /**
     * Returns the cross-sectional area of the lattice element.
     * @return Cross-section area.
     */
    virtual double giveArea() = 0;

    /**
     * Returns the radius of the lattice element.
     * @return radius.
     */
    virtual double giveRadius(){return 0;}

    /**
     * Returns the element length
     * @return Element length.
     */
    virtual double giveLength() = 0;

    /**
     * Returns the pressure
     * @return pressure.
     */
    virtual double givePressure() { return 0; }

    /**
     * Returns the old pressure
     * @return old pressure.
     */
    virtual double giveOldPressure() { return 0; }

    /**
     * Returns the mass
     * @return mass.
     */
    virtual double giveMass() { return 0; }
    
    /**
     * Returns the coupling flag
     * @return couplingFlag.
     */
    virtual int giveCouplingFlag(){return 0;}

    /**
     * Returns the coupling numbers
     * @return couplingNumbers.
     */
    virtual void giveCouplingNumbers(IntArray &numbers){return;}

    /**
     * Returns the crack lengths
     * @return crackLengths.
     */
    virtual void giveCrackLengths(FloatArray &lengths){return;}

    /**
     * Returns the crack widths
     * @return crackWidths.
     */
    virtual void giveCrackWidths(FloatArray &widths) {return;}


    /**
     * Gives the GP coordinates
     */
    virtual void  giveGpCoordinates(FloatArray &coords) { return; }

    /**
     * Returns the real volume
     */
    virtual double  computeLatticeVolumeAroundIt() {return 0;}


};
} // end namespace oofem
#endif
